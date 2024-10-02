#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>

#include "../common/version.h"
#include "../common/types/satc_data.h"


struct Params {

	inline const static std::string default_output_ext = ".satc";
	inline const static std::string default_sample_name_to_id_mapping_ext = ".sn2id.txt";

	std::string input;
	std::string output;
	std::string sample_name_to_id_mapping;
	uint32_t anchor_len = 0;  //required if the dump is in splash format
	uint32_t gap_len = 0;
	uint64_t max_cnt = std::numeric_limits<uint16_t>::max(); //I think 2 bytes is a reasonable maximum
	std::string ordering = "SBATC";
	uint32_t num_bytes_sample_id = 2;

	void Print(std::ostream& oss) const {
		oss << "input                : " << input << "\n";
		oss << "output               : " << output << "\n";
		oss << "output mapping       : " << sample_name_to_id_mapping << "\n";
		oss << "anchor len           : " << anchor_len << "\n";
		oss << "gap len              : " << gap_len << "\n";
		oss << "max cnt              : " << max_cnt << "\n";
		oss << "ordering             : " << ordering << "\n";
		oss << "num_bytes_sample_id  : " << num_bytes_sample_id << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_undump\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options]\n";

		Params defaults;
		std::cerr
			<< "Options:\n"
			<< "    -i <path> [REQUIRED] - path to input file (that is dumped satc)" << "\n"
			<< "    -o <path>            - path to output SATC file (default: input path with added extension " << default_output_ext << ")" << "\n"
			<< "    -m <path>            - path to output sample name to id mapping (default: input path with added extension " << default_sample_name_to_id_mapping_ext << ")" << "\n"
			<< "    -c <int>             - max value of count to be stored in satc file, if any input value is larger it will be stored as this value (default: " << defaults.max_cnt <<")\n"
			<< "    --anchor_len <int>   - anchor_len, required if the dump is in splash format, ignored otherwise\n"
			<< "    --gap_len <int>      - gap_len, will be stored in satc file (default: " << defaults.gap_len << ")\n"
			<< "    -s <int>             - num bytes for sample id, number of unique samples in the input file must be <= pow(2, 8 * s) (default: " << defaults.num_bytes_sample_id << ")\n"
			<< "    -r <string>          - ordering, one of { ATSBC, SBATC, TASBC } - if input file is sorted in some order setting this appropriately may reduce output satc file size (default: " << defaults.ordering << ")\n"
			;
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

		if (param == "-i" && i + 1 < argc)
			res.input = argv[++i];
		else if (param == "-o" && i + 1 < argc)
			res.output = argv[++i];
		else if (param == "-m" && i + 1 < argc)
			res.sample_name_to_id_mapping = argv[++i];
		else if (param == "--anchor_len" && i + 1 < argc)
			res.anchor_len = std::stoul(argv[++i]);
		else if (param == "-c" && i + 1 < argc)
			res.max_cnt = std::stoul(argv[++i]);
		else if (param == "--gap_len" && i + 1 < argc)
			res.gap_len = std::stoul(argv[++i]);
		else if (param == "-s" && i + 1 < argc)
			res.num_bytes_sample_id = std::stoul(argv[++i]);
		else if (param == "-r" && i + 1 < argc)
			res.ordering = argv[++i];
		else {
			std::cerr << "Error: unknown option " << param << "\n";
			exit(1);
		}
	}

	if (res.input.empty()) {
		std::cerr << "Error: input missing (-i)\n";
		exit(1);
	}

	if (res.output.empty())
		res.output = res.input + Params::default_output_ext;

	if (res.sample_name_to_id_mapping.empty())
		res.sample_name_to_id_mapping = res.input + Params::default_sample_name_to_id_mapping_ext;

	return res;
}


void GuessConfig(const std::string& path, RecFmt& fmt, uint32_t& anchor_len, uint32_t& target_len, uint32_t& barcode_len) {
	
	std::string line;
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open file " << path << "\n";
		exit(1);
	}
	if (!std::getline(in, line)) {
		std::cerr << "Error: cannot read first line from file" << path << "\n";
		exit(1);
	}

	//in SATC format we have <sample> [barcode] <anchor> <target> <count> -> so 4 or 5 parts
	//in SPLASH format we have <count> <anchor+target> <sample[_barcode]> -> so 3 parts

	std::vector<std::string> parts;
	std::string tmp;
	std::istringstream iss(line);
	while (iss >> tmp)
		parts.emplace_back(tmp);

	auto is_dna_seq = [](const std::string& seq) -> bool
		{
			return seq.find_first_not_of("ACGT") == std::string::npos;
		};

	std::string empty_string;

	std::string sample, barcode, anchor, target, count;

	if (parts.size() == 4 || parts.size() == 5) {
		fmt = RecFmt::SATC;
		sample = parts[0];

		//if we have barcode
		if (parts.size() == 5) {
			barcode = parts[1];
			anchor = parts[2];
			target = parts[3];
			count = parts[4];
		}
		else {
			anchor = parts[1];
			target = parts[2];
			count = parts[3];
		}
		

	}
	else if (parts.size() == 3) {
		fmt = RecFmt::SPLASH;
		if (anchor_len == 0) {
			std::cerr << "Error: in case of SPLASH dump format anchor len must be provided\n";
			exit(1);
		}

		count = parts[0];
		std::string anchor_target = parts[1];
		std::string sample_barcode = parts[2];

		if (anchor_len > anchor_target.length()) {
			std::cerr << "Error: invalid anchor len provided, anchor + targets are of length: " << anchor_target.length() << "\n";
			exit(1);
		}

		anchor = anchor_target.substr(0, anchor_len);
		target = anchor_target.substr(anchor_len);

		auto p = sample_barcode.find_last_of('_');
		if (p != std::string::npos) {
			auto maybe_barcode = sample_barcode.substr(p + 1);
			if (maybe_barcode.length() > 6 && is_dna_seq(maybe_barcode)) //because maybe we have a sample name ABC_1, meaning 1 is just a part of regular name, I will assume that if after _ we have at least 6 characters and only ACGT this is barcode
				barcode = maybe_barcode;
		}
	}
	else {
		std::cerr << "Error: cannot determine format from line " << line << "\n";
		exit(1);
	}


	barcode_len = barcode.length();
	anchor_len = anchor.length();
	target_len = target.length();

	if (!is_dna_seq(barcode)) {
		std::cerr << "Error: invalid barcode: " << barcode << "\n";
		exit(1);
	}

	if (!is_dna_seq(anchor)) {
		std::cerr << "Error: invalid anchor: " << barcode << "\n";
		exit(1);
	}

	if (!is_dna_seq(target)) {
		std::cerr << "Error: invalid target: " << barcode << "\n";
		exit(1);
	}
}

class SampleNameToIdMapping {
	std::unordered_map<std::string, uint32_t> map;
public:
	uint32_t get(const std::string& sample_name) {
		return map.emplace(sample_name, static_cast<uint32_t>(map.size())).first->second;
	}

	std::vector<std::string> GetSortedByIdMapping() const {
		std::vector<std::string> res(map.size());
		for (const auto& [name, id] : map)
			res[id] = name;
		return res;
	}
};

class SATC_dump_record_walker {
	std::ifstream in;
	uint32_t anchor_len;

	RecFmt fmt{};
	uint32_t target_len{};
	uint32_t barcode_len{};

	std::string line;
	std::string tmp;

	SampleNameToIdMapping sample_name_to_id_mapping;
	//mkokot_TODO: consider not using istringstream...

	void Extract_SATC_fmt(Record& rec) {
		std::istringstream iss(line);
		if (!(iss >> tmp)){
			std::cerr << "Error: cannot read sample from line: " << line << "\n";
			exit(1);
		}
		rec.sample_id = sample_name_to_id_mapping.get(tmp);
		if (barcode_len) {
			if (!(iss >> tmp)) {
				std::cerr << "Error: cannot read barcode from line: " << line << "\n";
				exit(1);
			}
			if (tmp.length() != barcode_len) {
				std::cerr << "Error: inconsistent barcode lengths, detected in line " << line << "\n";
				exit(1);
			}
			rec.barcode = str_kmer_to_uint64_t(tmp);
		}

		if (!(iss >> tmp)) {
			std::cerr << "Error: cannot read anchor from line: " << line << "\n";
			exit(1);
		}
		if (tmp.length() != anchor_len) {
			std::cerr << "Error: inconsistent anchor lengths, detected in line " << line << "\n";
			exit(1);
		}
		rec.anchor = str_kmer_to_uint64_t(tmp);

		if (!(iss >> tmp)) {
			std::cerr << "Error: cannot read target from line: " << line << "\n";
			exit(1);
		}
		if (tmp.length() != target_len) {
			std::cerr << "Error: inconsistent target lengths, detected in line " << line << "\n";
			exit(1);
		}
		rec.target = str_kmer_to_uint64_t(tmp);
		if (!(iss >> rec.count)) {
			std::cerr << "Error: cannot read count from line: " << line << "\n";
			exit(1);
		}
	}

	void Extract_SPLASH_fmt(Record& rec) {
		std::istringstream iss(line);
		if (!(iss >> rec.count)) {
			std::cerr << "Error: cannot read count from line: " << line << "\n";
			exit(1);
		}

		if (!(iss >> tmp)) {
			std::cerr << "Error: cannot read anchor+target from line: " << line << "\n";
			exit(1);
		}

		if (tmp.length() != anchor_len + target_len) {
			std::cerr << "Error: inconsistent anchor+target lengths, detected in line " << line << "\n";
			exit(1);
		}
		rec.anchor = str_kmer_to_uint64_t(tmp.substr(0, anchor_len));
		rec.target = str_kmer_to_uint64_t(tmp.substr(anchor_len));

		if (!(iss >> tmp)) {
			std::cerr << "Error: cannot read sample from line: " << line << "\n";
			exit(1);
		}

		if (barcode_len) {
			
			auto sample_name = tmp.substr(0, tmp.length() - barcode_len);
			auto barcode = tmp.substr(tmp.length() - barcode_len);
			rec.barcode = str_kmer_to_uint64_t(barcode);
			rec.sample_id = sample_name_to_id_mapping.get(sample_name);
		}
		else {
			rec.sample_id = sample_name_to_id_mapping.get(tmp);
		}
	}

public:
	SATC_dump_record_walker(
		const std::string& path,
		uint32_t anchor_len
		) :
		in(path),
		anchor_len(anchor_len)
	{
		if (!in)
		{
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		GuessConfig(path, fmt, this->anchor_len, target_len, barcode_len);
	}

	uint32_t AnchorLen() const
	{
		return anchor_len;
	}

	uint32_t TargetLen() const
	{
		return target_len;
	}

	uint32_t BarcodeLen() const
	{
		return barcode_len;
	}

	bool Next(Record& rec) {
		if (!std::getline(in, line))
			return false;
		if (fmt == RecFmt::SATC)
			Extract_SATC_fmt(rec);
		else if (fmt == RecFmt::SPLASH)
			Extract_SPLASH_fmt(rec);
		else
			assert(false);

		return true;
	}

	std::vector<std::string> GetSampleNameMapping() const {
		return sample_name_to_id_mapping.GetSortedByIdMapping();
	}
};

//mkokot_TODO: copied from fafq_filter/worker.h -> refactor this
constexpr uint8_t no_bytes(size_t x)
{
	uint8_t nb = 1;
	x /= 256;

	while (x)
	{
		++nb;
		x /= 256;
	}

	return nb;
}

void process(const Params& params) {
	buffered_binary_writer satc_writer(params.output);
	if (!satc_writer) {
		std::cerr << "Error: cannot open output file " << params.output << "\n";
		exit(1);
	}

	SATC_dump_record_walker walker(params.input, params.anchor_len);

	Header header;
	assert(walker.AnchorLen() <= 255);
	assert(walker.TargetLen()<= 255);
	assert(walker.BarcodeLen()<= 255);
	assert(params.gap_len <= 255);
	assert(params.num_bytes_sample_id <= 255);
	
	header.sample_id_size_bytes = params.num_bytes_sample_id;
	
	header.counter_size_bytes = no_bytes(params.max_cnt);
	header.gap_len_symbols = static_cast<uint8_t>(params.gap_len);

	header.anchor_len_symbols = static_cast<uint8_t>(walker.AnchorLen());
	header.target_len_symbols = static_cast<uint8_t>(walker.TargetLen());
	header.barcode_len_symbols = static_cast<uint8_t>(walker.BarcodeLen());

	header.anchor_size_bytes = (header.anchor_len_symbols + 3) / 4;
	header.target_size_bytes = (header.target_len_symbols + 3) / 4;
	header.barcode_size_bytes = (header.barcode_len_symbols + 3) / 4;

	header.ordering = Header::ordering_from_string(params.ordering);

	header.serialize(satc_writer);

	auto max_sample_id = (1ul << (8 * params.num_bytes_sample_id)) - 1;

	Record rec;
	while (walker.Next(rec)) {
		if (rec.count > params.max_cnt)
			rec.count = params.max_cnt;
		if (rec.sample_id > max_sample_id) {
			std::cerr << "Error: number of unique samples is too large for " << params.num_bytes_sample_id << " bytes used to store sample id, rerun with increased -s\n";
			exit(1);
		}
		rec.serialize(satc_writer, header);
	}

	std::ofstream sample_name_to_id(params.sample_name_to_id_mapping);
	if (!sample_name_to_id) {
		std::cerr << "Error: cannot open output file " << params.sample_name_to_id_mapping << "\n";
		exit(1);
	}

	auto mapping = walker.GetSampleNameMapping();

	for (size_t id = 0 ; id < mapping.size() ; ++id) {
		sample_name_to_id << mapping[id] << ' ' << id << "\n";
	}
}

int main(int argc, char**argv) {
	auto params = read_params(argc, argv);
	params.Print(std::cerr);

	process(params);

	return 0;
}