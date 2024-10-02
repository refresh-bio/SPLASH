#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <regex>
#include <cstring>
#include <cassert>
#include <sstream>
#include <limits>
#include <algorithm>
#include <cstdint>
#include  "../common/version.h"

//mkokot_TODO: this is copied from lookup_table
//quite ugly because similar class is in satc_data, but it was extended for compression
//so here is simplified variant
class simple_buffered_binary_writer {
	std::ofstream out;
	std::vector<uint8_t> buff;

	void flush() {
		out.write(reinterpret_cast<char*>(buff.data()), buff.size());
		buff.clear();
	}

	void assure_space(size_t size) {
		if (buff.size() + size > buff.capacity()) {
			flush();
			if (size > buff.capacity())
				buff.reserve(size);
		}
	}

	void write(const std::vector<uint8_t>& vec) {
		assure_space(vec.size());
		for (auto x : vec)
			buff.push_back(x);
	}

public:
	simple_buffered_binary_writer(simple_buffered_binary_writer&&) = default;
	simple_buffered_binary_writer& operator=(simple_buffered_binary_writer&& rhs) noexcept {
		if (&rhs == this)
			return *this;
		close();

		out = std::move(rhs.out);
		buff = std::move(rhs.buff);

		return *this;
	}
	simple_buffered_binary_writer(const std::string& path, size_t buff_size = 1ull << 25) {
		out.rdbuf()->pubsetbuf(nullptr, 0);
		out.open(path, std::ios::binary);
		buff.reserve(buff_size);
	}
	operator bool() const {
		return out.operator bool();
	}
	void write(const uint8_t* ptr, size_t size) {
		assure_space(size);
		buff.insert(buff.end(), ptr, ptr + size);
	}

	void write(const char* ptr, size_t size) {
		write(reinterpret_cast<const uint8_t*>(ptr), size);
	}

	void write(const char* ptr) {
		write(ptr, strlen(ptr));
	}

	void write(const std::string& str) {
		write(str.c_str(), str.size());
	}

	void write(char c)
	{
		assure_space(1);
		buff.push_back(static_cast<uint8_t>(c));
	}

	~simple_buffered_binary_writer() {
		if (!buff.empty())
			flush();
	}
	void close() {
		if (!buff.empty())
			flush();
		out.close();
	}
};

//mkokot_TODO: based on implementation in lookup tables, consider refactor
//mkokot_TODO: consider resignation from streams...
class SeqReaderExtendors {
	std::vector<std::string> col_names;
	std::string anchor;
	std::vector<std::string> targets;
	size_t cur_target_id{};
	bool finished = false;

	std::vector<size_t> target_cols_ids;
	std::ifstream in;
	uint32_t n_most_freq_targets; //only for extendors, skipped for compactors
	uint32_t cur_writen_targets;
	std::string line;

	bool is_col_with_target(const std::string& col_name) {
		const static std::regex pattern("most_freq_target_([1-9][0-9]*)");
		return std::regex_match(col_name, pattern);
	}

	bool init_extendor() {
		cur_target_id = 0;
		cur_writen_targets = 0;

		if (!(std::getline(in, line)) || line == "")
			return false;

		std::istringstream iss(line);

		if (!(iss >> anchor)) {
			std::cerr << "Error: cannot read anchor from line " << line;
			exit(1);
		}
		size_t id = 1;
		size_t i = 0;
		std::string val;
		while (iss >> val && i < target_cols_ids.size()) {
			if (id == target_cols_ids[i]) {
				targets[i++] = val;
			}
			++id;
		}

		if (i != target_cols_ids.size()) {
			std::cerr << "Error: cannot read all targets from line " << line << "\n";
			exit(1);
		}

		return true;
	}

public:
	SeqReaderExtendors(const std::string& path, uint32_t n_most_freq_targets) :
		in(path),
		n_most_freq_targets(n_most_freq_targets) {

		assert(n_most_freq_targets);
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		if (!std::getline(in, line)) {
			std::cerr << "Error: cannot read header from file " << path << "\n";
			exit(1);
		}
		std::istringstream iss(line);

		std::string name;
		while (iss >> name)
			col_names.emplace_back(std::move(name));

		if (col_names[0] != "anchor") {
			std::cerr << "Error: first columns should be \"anchor\" but is " << col_names[0] << "\n";
			exit(1);
		}

		for (size_t i = 0; i < col_names.size(); ++i) {
			if (is_col_with_target(col_names[i]))
				target_cols_ids.push_back(i);
		}

		if (target_cols_ids.empty()) {
			std::cerr << "Error: this file does not contain targets (non of columns start with pattern \"most_freq_target_\"\n";
			exit(1);
		}

		targets.resize(target_cols_ids.size());

		finished = !init_extendor();
	}

	bool NextSeq(simple_buffered_binary_writer& writer) {
		while (true) {

			if (finished)
				return false;

			if (cur_target_id < targets.size() && cur_writen_targets < n_most_freq_targets) {
				const std::string& target = targets[cur_target_id++];
				if (target == "-")
					continue;

				writer.write(">\n", 2);
				writer.write(anchor);
				writer.write(target);
				writer.write("\n", 1);
				++cur_writen_targets;
				return true;
			}

			finished = !init_extendor();
		}
	}
};

class SeqReaderCompactors {
	std::ifstream in;

	std::string header;
	std::string current_line;
	std::string compactor_seq;

	size_t compactor_pos = std::numeric_limits<size_t>::max();

public:
	SeqReaderCompactors(const std::string& path) :in(path) {
		if (!in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}

		if (!std::getline(in, header)) {
			std::cerr << "Error: cannot read header from file " << path << "\n";
			exit(1);
		}

		std::istringstream iss(header);

		size_t i{};
		std::string header_part;
		while (iss >> header_part) {
			if (header_part == "compactor") {
				compactor_pos = i;
				break;
			}
			++i;
		}

		if (compactor_pos == std::numeric_limits<size_t>::max()) {
			std::cerr << "Error: cannot find column with name \"compactor\"\n";
			exit(1);
		}
	}
	bool NextSeq(simple_buffered_binary_writer& writer) {

		if (!std::getline(in, current_line))
			return false;

		if (current_line == "")
			return false;

		std::istringstream iss(current_line);

		size_t i{};
		while (iss >> compactor_seq) {
			if (i == compactor_pos) {
				writer.write(">\n", 2);
				writer.write(compactor_seq);
				writer.write("\n", 1);
				return true;
			}
			++i;
		}

		std::cerr << "Error: cannot read compactor from line " << current_line << "\n";
		exit(1);
	}

};

struct Params
{
	uint32_t n_most_freq_targets = std::numeric_limits<uint32_t>::max();

	std::string input;
	std::string output;
	void Print(std::ostream& oss) const
	{
		oss << "n_most_freq_targets            : " << n_most_freq_targets << "\n";
		oss << "input                          : " << input << "\n";
		oss << "output                         : " << output << "\n";
	}
	static void Usage(char* prog_name) {
		std::cerr << "tsv_to_fasta\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <input> <output>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input>      - path to splash result in tsv format containing anchors and targets\n"
			<< "    <output>     - path to output fasta file with extendors\n";
		std::cerr
			<< "Options:\n"
			<< "    --n_most_freq_targets <int> - how many (at most) most freq targets to take to build extendor for each anchor\n";
	}
};

Params read_params(int argc, char** argv)
{
	Params res;
	if (argc == 1) {
		Params::Usage(argv[0]);
		exit(0);
	}
	int i = 1;
	for (; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;

		std::string param = argv[i];
		if (param == "--n_most_freq_targets") {
			std::string tmp = argv[++i];
			res.n_most_freq_targets = std::stoul(tmp);
			if (res.n_most_freq_targets == 0) {
				std::cerr << "Error: --n_most_freq_targets must be >0\n";
				exit(1);
			}
		}
	}
	
	if (i >= argc) {
		std::cerr << "Error: input path must be specified\n";
		exit(1);
	}
	res.input = argv[i++];
	
	if (i >= argc) {
		std::cerr << "Error: output path must be specified\n";
		exit(1);
	}
	res.output = argv[i++];

	return res;
}

enum class InputFormat {extendors, copactors, undetermined};
InputFormat determine_input_type(const std::string& path)
{
	std::ifstream in(path);
	if (!in)
	{
		std::cerr << "Error: cannot open file " << path << "\n";
		exit(1);
	}
	std::string line;
	if (!std::getline(in, line)) {
		std::cerr << "Error: cannot read header from file " << path << "\n";
		exit(1);
	}
	std::istringstream iss(line);

	std::vector<std::string> col_names;
	std::string name;
	while (iss >> name)
		col_names.emplace_back(std::move(name));

	if (col_names.empty())
	{
		std::cerr << "Error: cannot extract column names from file " << path << "\n";
		exit(1);
	}
	if (col_names[0] != "anchor")
	{
		std::cerr << "Error: unsupported format, first column name should be 'anchor' in file: " << path << "\n";
		exit(1);
	}
	for (const auto& col_name : col_names)
	{
		if (col_name == "compactor")
			return InputFormat::copactors;
		if (col_name.find("most_freq_target_") == 0)
			return InputFormat::extendors;
	}

	return InputFormat::undetermined;
}

int main(int argc, char** argv) {
	std::cerr << "Welcome to tsv_to_fasta - utility to convert SPLASH tsv output file to FASTA with extendors/compactors\n";
	auto params = read_params(argc, argv);

	params.Print(std::cerr);
	auto input_format = determine_input_type(params.input);
	if (input_format == InputFormat::undetermined)
	{
		std::cerr << "Error: cannot determine input file format for: " << params.input << "\n";
		return 1;
	}
	if (input_format == InputFormat::extendors)
	{
		simple_buffered_binary_writer writer(params.output);
		SeqReaderExtendors reader(params.input, params.n_most_freq_targets);
		while (reader.NextSeq(writer))
			;
	}
	else if (input_format == InputFormat::copactors)
	{
		simple_buffered_binary_writer writer(params.output);
		SeqReaderCompactors reader(params.input);
		while (reader.NextSeq(writer))
			;
	}
}
