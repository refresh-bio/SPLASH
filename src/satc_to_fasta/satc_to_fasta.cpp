#define _CRT_SECURE_NO_WARNINGS

// (option for the script: samples in this list only if > L of the anchors are detected in the sample) - to moze powinno byc w SATC filter?
// https://canlifeexist.slack.com/archives/C044H6V26DB/p1716649451039259
// no i moze jakies jeszcze tryby
// sprawdzamy czy w mappingach nie ma duplikatow
// podczas ladowania satcy sprawdzamy czy dane sa juz posortowane tak jak tego potrzebujemy
// potem mamy posortowane
// mamy tez liste anchorow

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <unordered_map>

#include "../common/version.h"
#include "../common/accepted_anchors.h"

struct Params {

	struct InputDesc
	{
		std::string satc_path;
		std::string sample_name_to_id_mapping_path;
	};
	std::vector<InputDesc> input;
	std::string output;

	std::string anchor_list;

	//uint32_t t = (std::min)(8u, std::thread::hardware_concurrency()); //mkokot_TODO: parallelize this software?
	uint32_t n = 0; //number of targets to be printed per anchor, if there is more in the SATC file, only the n most abundant will be printed, if there is less (including zero if list of anchors provided) the missing one are replaced with Ns
	bool l = false;
	bool s = false; //skip anchors (print only targets)

	void Print(std::ostream& oss) const
	{
		oss << "input         :\n";
		for (const auto& x : input)
			oss << "\t" << x.satc_path << "\t" << x.sample_name_to_id_mapping_path << "\n";
		oss << "output        : " << output << "\n";
		oss << "anchor_list   : " << anchor_list << "\n";
		oss << "n             : " << n << "\n";
		oss << "l             : " << std::boolalpha << l << "\n";
		oss << "s             : " << std::boolalpha << s << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_to_fasta\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options]\n";

		Params defaults;
		std::cerr
			<< "Options:\n"
			<< "    -i [<path_satc> <path_sample_names>|<path>]\n"
			"              [REQUIRED] - paths to input SATC file and sample name (or txt with list of paths, if -l used), may be used multiple times" << "\n"
			<< "    -o <path> [REQUIRED] - path to output FASTA file \n"
			<< "    -d <path>            - path to file with anchors dictionary (k-mer per line), if not provided will just print anchors from SATC, if provided it must be assured that all anchors in SATC files are also on this list, but if anchor from this list is not in SATC file for given sample it will be printed anyway with Ns as target" << "\n"
			//			<< "    -t <int>             - number of threads, default: " << defaults.t << "\n"
			<< "    -n <int>             - number of most freq targets to print for each anchor in each sample, if there is more in the SATC file, only the n  most abundant will be printed, if there is less (including zero if list of anchors provided) the missing one are replaced with Ns, usually should be set the same as for satc_filter, default: " << defaults.n << "\n"
			<< "    -l                   - if set <input> is treated as lists of file paths\n"
			<< "    -s                   - skip anchors (only targets are printed)\n";
	}
};

bool contains_flag(int argc, char** argv, const std::string& flag)
{
	for (int i = 1; i < argc; ++i)
		if (argv[i] == flag)
			return true;
	return false;
}

void read_input_list(const std::string& path, std::vector<Params::InputDesc>& res) {
	res.clear();
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open file " << path << "\n";
		exit(1);
	}
	std::string line;
	Params::InputDesc item;
	while (std::getline(in, line)) {
		if (line.empty())
			continue;

		std::istringstream iss(line);
		if (! (iss >> item.satc_path >> item.sample_name_to_id_mapping_path)) {
			std::cerr << "Error: cannot read <path_satc> <path_sample_names> from line " << line << " in " << path << "\n";
			exit(1);
		}

		if ((iss >> line) && !line.empty()) {
			std::cerr << "Error: unexpected data in line " << line << " in " << path << "\n";
			std::cerr << "(line should contain only a pair of whitespace separated pair of paths, whitespaces inside paths are not allowed\n";
			exit(1);
		}
		res.emplace_back(std::move(item));
	}
}

Params read_params(int argc, char** argv)
{
	Params res;

	if (argc == 1) {
		Params::Usage(argv[0]);
		exit(0);
	}

	//must know if first to know how to tread -i (either pair of paths or path to a file with list of pair of paths)
	res.l = contains_flag(argc, argv, "-l");

	int i = 1;
	for (; i < argc; ++i) {
		if (argv[i][0] != '-')
			break;

		std::string param = argv[i];

		if (param == "-i" && i + 1 < argc) {

			if (res.l)
				read_input_list(argv[++i], res.input);
			else {
				res.input.emplace_back();
				res.input.back().satc_path = argv[++i];
				if (i + 1 < argc) {
					res.input.back().sample_name_to_id_mapping_path = argv[++i];
					if (res.input.back().sample_name_to_id_mapping_path.empty()) {
						std::cerr << "Error: <path_sample_names> is required for -i when -l not used\n";
						exit(1);
					}
					if (res.input.back().sample_name_to_id_mapping_path[0] == '-') {
						std::cerr << "Error: <path_sample_names> should be a path but " << res.input.back().sample_name_to_id_mapping_path << " found\n";
						exit(1);
					}
				}
				else {
					std::cerr << "Error: <path_sample_names> is required for -i when -l not used\n";
					exit(1);
				}
			}
		}
		else if (param == "-o" && i + 1 < argc)
			res.output = argv[++i];
		else if (param == "-d" && i + 1 < argc)
			res.anchor_list = argv[++i];
		else if (param == "-n" && i + 1 < argc)
			res.n = std::stoul(argv[++i]);
		else if (param == "-s")
			res.s = true;

		//		else if (param == "-t" && i + 1 < argc)
		//			res.t = std::stoul(argv[++i]);

		else if (param == "-l") //-l was captured earlier but we need this to prevent treating this as unknown parameter
			;
		else {
			std::cerr << "Error: unknown option " << param << "\n";
			exit(1);
		}
	}

	//if (res.t == 0) {
	//	res.t = Params{}.t; //get default
	//}

	if (res.input.empty()) {
		std::cerr << "Error: input missing (-i)\n";
		exit(1);
	}

	if (res.output.empty()) {
		std::cerr << "Error: output missing (-o)\n";
		exit(1);
	}

	return res;
}

void verify_sample_names_unique(const std::vector<Params::InputDesc>& input) {
	using SampleNameType = std::string;
	using PathType = std::string;
	std::unordered_map<SampleNameType, PathType> map; //key - sample_name, value - 
	for (const auto& x : input) {
		const auto& mapping_path = x.sample_name_to_id_mapping_path;
		std::ifstream in(mapping_path);
		if (!in) {
			std::cerr << "Error: cannot open file " << mapping_path << "\n";
			exit(1);
		}
		SampleNameType sample_name;
		int id;
		while (in >> sample_name >> id) {
			auto it = map.find(sample_name);
			if (it != map.end()) {
				std::cerr << "Error: sample " << sample_name << " is in " << mapping_path << ", but was also in " << it->second << " - this is not allowed\n";
				exit(1);
			}
			else
				map.emplace_hint(it, std::make_pair(sample_name, mapping_path));
		}
	}
}

class StoreSeq {
	std::string cur_seq;
	uint32_t max_line_len = 60;
	//uint32_t max_line_len = 54; //debug purposes

	std::ofstream& out;
public:
	StoreSeq(std::ofstream& out) :
		out(out)
	{

	}
	void Add(const std::string& seq) {

		size_t i = 0;
		while (cur_seq.length() + seq.length() - i > max_line_len) {
			auto from_seq = max_line_len - cur_seq.length();
			cur_seq.insert(cur_seq.end(), seq.begin() + i, seq.begin() + i + from_seq);
			cur_seq.push_back('\n');
			out.write(cur_seq.c_str(), cur_seq.size());
			cur_seq.clear();
			i += from_seq;
		}

		cur_seq.insert(cur_seq.end(), seq.begin() + i, seq.end());
	}
	~StoreSeq() {
		if (cur_seq.length()) {
			cur_seq.push_back('\n');
			out.write(cur_seq.c_str(), cur_seq.size());
		}
	}
};

std::vector<std::string> load_mapping(const std::string& path) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open file " << path << "\n";
		exit(1);
	}
	std::vector<std::string> res;
	std::string name;
	size_t id;
	while (in >> name >> id) {
		if (id >= res.size()) {
			res.resize(id + 1);
		}
		res[id] = name;
	}
	return res;
}

struct SampleAnchorTargetCount {
	uint32_t sample_id;
	uint64_t anchor;
	uint64_t target;
	uint32_t count;
	SampleAnchorTargetCount(
		uint32_t sample_id,
		uint64_t anchor,
		uint64_t target,
		uint32_t count) :
		sample_id(sample_id),
		anchor(anchor),
		target(target),
		count(count)
	{}
	bool operator<(const SampleAnchorTargetCount& rhs) const {
		//count from largest
		return std::make_tuple(sample_id, anchor, rhs.count, target) <
			std::make_tuple(rhs.sample_id, rhs.anchor, count, rhs.target);
	}
};

struct SampleBarcodeAnchorTargetCount {
	uint32_t sample_id;
	uint64_t barcode;
	uint64_t anchor;
	uint64_t target;
	uint32_t count;
	SampleBarcodeAnchorTargetCount(
		uint32_t sample_id,
		uint64_t barcode,
		uint64_t anchor,
		uint64_t target,
		uint32_t count) :
		sample_id(sample_id),
		barcode(barcode),
		anchor(anchor),
		target(target),
		count(count)
	{}
	bool operator<(const SampleBarcodeAnchorTargetCount& rhs) const {
		//count from largest
		return std::make_tuple(sample_id, barcode, anchor, rhs.count, target) <
			std::make_tuple(rhs.sample_id, rhs.barcode, rhs.anchor, count, rhs.target);
	}
};

//mkokot_TODO: this code is very similar to this in a function calling it (below the call)
//so this should be refactored...
void process_single_10X(
	const Params::InputDesc& input,
	uint32_t anchor_len,
	buffered_binary_reader& satc_reader,
	Header& header,
	const std::vector<std::string>& mapping,
	const std::vector<uint64_t>& sorted_anchors,
	std::ofstream& out,
	bool skip_anchors,
	uint32_t n //mkokot_TODO: use this parameter
	) {
	std::vector<SampleBarcodeAnchorTargetCount> data;

	Record rec;
	bool is_already_sorted = true;
	if (rec.load(satc_reader, header)) {
		data.emplace_back(rec.sample_id, rec.barcode, rec.anchor, rec.target, rec.count);
		while (rec.load(satc_reader, header)) {
			data.emplace_back(rec.sample_id, rec.barcode, rec.anchor, rec.target, rec.count);
			is_already_sorted &= data[data.size() - 2] < data.back();
		}
	}
	if (!is_already_sorted)
		std::sort(data.begin(), data.end());

	std::string Ns(header.target_len_symbols, 'N');

	std::string anchor_str;

	//go through samples
	for (uint64_t i = 0; i < data.size();) {

		auto sample_id = data[i].sample_id;
		auto barcode = data[i].barcode;
		auto sample_name = mapping[sample_id];

		if (sample_name.empty()) {
			std::cerr << "Error: not sample name for sample id " << sample_id << "\n";
			exit(1);
		}

		auto fasta_header = ">" + sample_name + " " + kmer_to_string(barcode, header.barcode_len_symbols) + "\n";
		out.write(fasta_header.c_str(), fasta_header.length());

		StoreSeq store_seq(out);

		size_t sorted_anchors_idx = 0;

		//go through anchors for given sample
		while (i < data.size() && data[i].sample_id == sample_id && data[i].barcode == barcode) {

			auto anchor = data[i].anchor;

			anchor_str = kmer_to_string(anchor, anchor_len);

			if (!sorted_anchors.empty()) {
				//at some point we must find anchor from data in accepted anchors
				while (sorted_anchors_idx < sorted_anchors.size() && sorted_anchors[sorted_anchors_idx] != anchor) {
					if (sorted_anchors[sorted_anchors_idx] > anchor) {
						std::cerr << "Error: SATC file " << input.satc_path << " contains anchor (" << anchor_str << ") that is not on provided anchors list - this is not supported\n";
						exit(1);
					}

					if (!skip_anchors)
						store_seq.Add(kmer_to_string(sorted_anchors[sorted_anchors_idx], anchor_len));

					for (uint32_t j = 0; j < n ; ++j)
						store_seq.Add(Ns);

					++sorted_anchors_idx;
				}
				if (sorted_anchors_idx == sorted_anchors.size()) {
					std::cerr << "Error: SATC file " << input.satc_path << " contains anchor (" << anchor_str << ") that is not on provided anchors list - this is not supported\n";
					exit(1);
				}
				++sorted_anchors_idx;
			}

			//go through targets
			if (!skip_anchors)
				store_seq.Add(anchor_str);

			std::string target_str;
			uint32_t num_printed_targets = 0;
			while (i < data.size() && data[i].sample_id == sample_id && data[i].barcode == barcode && anchor == data[i].anchor) {

				if (num_printed_targets < n) {
					target_str = kmer_to_string(data[i].target, header.target_len_symbols);
					store_seq.Add(target_str);
					++num_printed_targets;
				}

				++i;
			}

			for ( ;num_printed_targets < n; ++num_printed_targets)
				store_seq.Add(Ns);
		}
		//remaining anchors that were not found in the data
		//this will have 0 iterations if sorted_anchors is empty
		for (; sorted_anchors_idx < sorted_anchors.size(); ++sorted_anchors_idx) {
			if (!skip_anchors)
				store_seq.Add(kmer_to_string(sorted_anchors[sorted_anchors_idx], anchor_len));

			for (uint32_t j = 0; j < n ; ++j)
				store_seq.Add(Ns);
		}
	}
}

void process_single(
	const Params::InputDesc& input,
	uint32_t& anchor_len,
	const std::vector<uint64_t>& sorted_anchors, //may be empty
	std::ofstream& out,
	bool skip_anchors,
	uint32_t n) {

	Header header;
	buffered_binary_reader satc_reader(input.satc_path);
	if (!satc_reader) {
		std::cerr << "Error: cannot open SATC file " << input.satc_path << "\n";
		exit(1);
	}

	auto mapping = load_mapping(input.sample_name_to_id_mapping_path);

	header.load(satc_reader);

	if (anchor_len == 0)
		anchor_len = header.anchor_len_symbols;

	if (header.anchor_len_symbols != anchor_len) {
		std::cerr << "Error: inconsistent anchor lengths\n";
		exit(1);
	}

	if (header.barcode_len_symbols != 0) {
		process_single_10X(
			input,
			anchor_len,
			satc_reader,
			header,
			mapping,
			sorted_anchors,
			out,
			skip_anchors,
			n);
		return;
	}

	std::vector<SampleAnchorTargetCount> data;

	Record rec;
	bool is_already_sorted = true;
	if (rec.load(satc_reader, header)) {
		data.emplace_back(rec.sample_id, rec.anchor, rec.target, rec.count);
		while (rec.load(satc_reader, header)) {
			data.emplace_back(rec.sample_id, rec.anchor, rec.target, rec.count);
			is_already_sorted &= data[data.size() - 2] < data.back();
		}
	}
	if (!is_already_sorted)
		std::sort(data.begin(), data.end());

	std::string Ns(header.target_len_symbols, 'N');

	std::string anchor_str;

	//go through samples
	for (uint64_t i = 0; i < data.size();) {

		auto sample_id = data[i].sample_id;

		auto sample_name = mapping[sample_id];

		if (sample_name.empty()) {
			std::cerr << "Error: not sample name for sample id " << sample_id << "\n";
			exit(1);
		}

		auto fasta_header = ">" + sample_name + "\n";
		out.write(fasta_header.c_str(), fasta_header.length());

		StoreSeq store_seq(out);

		size_t sorted_anchors_idx = 0;

		//go through anchors for given sample
		while (i < data.size() && data[i].sample_id == sample_id) {

			auto anchor = data[i].anchor;

			anchor_str = kmer_to_string(anchor, anchor_len);

			if (!sorted_anchors.empty()) {
				//at some point we must find anchor from data in accepted anchors
				while (sorted_anchors_idx < sorted_anchors.size() && sorted_anchors[sorted_anchors_idx] != anchor) {
					if (sorted_anchors[sorted_anchors_idx] > anchor) {
						std::cerr << "Error: SATC file " << input.satc_path << " contains anchor (" << anchor_str << ") that is not on provided anchors list - this is not supported\n";
						exit(1);
					}

					if (!skip_anchors)
						store_seq.Add(kmer_to_string(sorted_anchors[sorted_anchors_idx], anchor_len));

					for (uint32_t j = 0; j < n ; ++j)
						store_seq.Add(Ns);

					++sorted_anchors_idx;
				}
				if (sorted_anchors_idx == sorted_anchors.size()) {
					std::cerr << "Error: SATC file " << input.satc_path << " contains anchor (" << anchor_str << ") that is not on provided anchors list - this is not supported\n";
					exit(1);
				}
				++sorted_anchors_idx;
			}

			//go through targets
			if (!skip_anchors)
				store_seq.Add(anchor_str);

			std::string target_str;
			uint32_t num_printed_targets = 0;
			while (i < data.size() && data[i].sample_id == sample_id && anchor == data[i].anchor) {

				if (num_printed_targets < n) {
					target_str = kmer_to_string(data[i].target, header.target_len_symbols);
					store_seq.Add(target_str);
					++num_printed_targets;
				}

				++i;
			}

			for ( ;num_printed_targets < n; ++num_printed_targets)
				store_seq.Add(Ns);
		}
		//remaining anchors that were not found in the data
		//this will have 0 iterations if sorted_anchors is empty
		for (; sorted_anchors_idx < sorted_anchors.size(); ++sorted_anchors_idx) {
			if (!skip_anchors)
				store_seq.Add(kmer_to_string(sorted_anchors[sorted_anchors_idx], anchor_len));

			for (uint32_t j = 0; j < n ; ++j)
				store_seq.Add(Ns);
		}
	}
}

void process(const Params& params) {

	verify_sample_names_unique(params.input); // mkokot_TODO: in fact this is to be considered, maybe we should allow sample names to be distributed to many files, but it have some implications that must be considered

	//now we could process each file in a separate thread, as we have distinct sample_names, but we would need to synchronize the output
	//via some queue
	//so lets start with a sequential version

	std::ofstream out(params.output, std::ios::binary);
	if (!out) {
		std::cerr << "Error: cannot open output file " << params.output << "\n";
		exit(1);
	}

	std::vector<uint64_t> sorted_accepted_anchors;

	uint32_t anchor_len = 0;

	std::vector<uint64_t> sorted_anchors;

	if (!params.anchor_list.empty()) {
		std::ifstream in(params.anchor_list);
		if (!in) {
			std::cerr << "Error: cannot open file " << params.anchor_list << "\n";
			exit(1);
		}

		ReadAnchorsFromPlainOrDSV reader(in, params.anchor_list, [&sorted_anchors](uint64_t anchor)
			{
				sorted_anchors.emplace_back(anchor);
			});

		anchor_len = reader.GetAnchorLen();

		std::sort(sorted_anchors.begin(), sorted_anchors.end());
	}



	for (const auto& x : params.input)
		process_single(x, anchor_len, sorted_anchors, out, params.s, params.n);
}

int main(int argc, char** argv) {
	auto params = read_params(argc, argv);
	params.Print(std::cerr);

	process(params);

	return 0;
}
