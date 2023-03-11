#include "../common/satc_data.h"
#include "../common/version.h"
#include "../common/accepted_anchors.h"
#include "../common/murmur64.h"
#include <string>
#include <iostream>
#include <vector>
#include <cstdint>
#include <unordered_set>

struct Params {
	std::string input;
	std::string output;
	std::string anchor_list_path;
	std::string sample_names;
	uint32_t n_bins = 0;
	bool separately = false;

	RecFmt format = RecFmt::SATC;
	void Print(std::ostream& oss) const
	{
		oss << "input         : " << input << "\n";
		oss << "output        : " << output << "\n";
		oss << "anchor_list   : " << anchor_list_path << "\n";
		oss << "sample_names  : " << sample_names << "\n";
		oss << "n_bins        : " << n_bins << "\n";
		oss << "separately    : " << std::boolalpha << separately << "\n";
		oss << "format        : " << RecFmtConv::to_string(format) << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_dump\n";
		NOMAD_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <input> <output>\n";
		std::cerr << "or\n";
		std::cerr << "\t" << prog_name << " --which-bin --n_bins <int> <anchor>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input> - path to binary file in satc format\n"
			<< "    <output> - output path\n";
		std::cerr
			<< "Options:\n"
			<< "    --anchor_list <path>  - path to text file containing anchors separated by whitespaces, only anchors from this file will be dumped\n"
			<< "    --sample_names <path> - path for decode sample id, each line should contain <sample_name> <sample_id>\n"
			<< "    --format <string>     - output format, available options: satc, nomad (default: satc)\n"
			<< "    --n_bins <int>        - if set to value different than 0 the input is interpreted as a list of bins (each bin in separate line, first list is bin_0, second line is bin_1, etc. (in case of ill-formed input results will be incorrect)\n"
			<< "    --separately          - if set with n_bins != 0 output param will be treated as suffix name and there will be output for each bin\n";
	}
};

bool which_bin_mode(int argc, char** argv) {
	uint32_t n_bins{};
	bool is_which_bin_mode = false;
	int i = 1;
	for (; i < argc; ++i)
	{
		if (argv[i][0] != '-')
			break;

		std::string param = argv[i];

		if (param == "--which-bin")
			is_which_bin_mode = true;
		else if (param == "--n_bins")
			n_bins = std::atol(argv[++i]);
	}
	if (!is_which_bin_mode)
		return false;

	if (i >= argc) {
		std::cerr << "Error: anchor sequence missing\n";
		exit(1);
	}
	std::string str_anchor = argv[i++];

	//std::cerr << "Anchor len: " << str_anchor.length() << "\n";
	auto anchor = str_kmer_to_uint64_t(str_anchor);
	uint64_t hash = MurMur64Hash{}(anchor);
	if (n_bins == 0) {
		std::cerr << "Warning: --n_bins was not specified, printng just a hash, take a mod <n_bins> yourself or rerun with --n_bins\n";
		std::cout << "MurMur64Hash: " << hash << "\n";
	}
	else {
		std::cout << "bin id: " << hash % n_bins << "\n";
	}

	return true;
}


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

		if (param == "--anchor_list")
			res.anchor_list_path = argv[++i];
		else if (param == "--sample_names")
			res.sample_names = argv[++i];
		else if (param == "--format")
			res.format = RecFmtConv::from_string(argv[++i]);
		else if (param == "--separately")
			res.separately = true;
		else if (param == "--n_bins")
			res.n_bins = std::atol(argv[++i]);
	}
	if (i >= argc) {
		std::cerr << "Error: input missing\n";
		exit(1);
	}
	res.input = argv[i++];

	if (i >= argc) {
		std::cerr << "Error: output missing\n";
		exit(1);
	}
	res.output = argv[i++];

	return res;
}

void process_single_bin_mode(const Params& params) {
	buffered_binary_reader in(params.input);
	if (!in) {
		std::cerr << "Error: cannot open file " << params.input << "\n";
		exit(1);
	}

	std::ofstream out(params.output);
	if (!out) {
		std::cerr << "Error: cannot open file " << params.output << "\n";
		exit(1);
	}

	Header header;
	header.load(in);
	header.print(std::cerr);

	AcceptedAnchors accepted_anchors(params.anchor_list_path);
	Record rec;
	SampleNameDecoder sample_name_decoder(params.sample_names);
	while (rec.load(in, header)) {
		if (accepted_anchors.IsAccepted(rec.anchor))
			rec.print(out, header, params.format, sample_name_decoder);
	}
}

std::vector<std::string> read_bins_paths(const std::string& path, uint32_t n_bins) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open file " << path;
		exit(1);
	}
	std::vector<std::string> bin_paths;
	std::string line;
	while (std::getline(in, line)) {
		if (line != "")
			bin_paths.emplace_back(std::move(line));
	}
	if (bin_paths.size() != n_bins) {
		std::cerr << "Error: the number of bin paths in " << path << " is " << bin_paths.size() << ", but the --n_bins is set to " << n_bins << "\n";
		exit(1);
	}
	return bin_paths;
}

std::vector<std::vector<uint64_t>> split_anchors(const std::string& path, uint32_t n_bins) {
	std::vector<std::vector<uint64_t>> res;
	if (path == "") {
		return res;
	}
	res.resize(n_bins);
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open file " << path << "\n";
		exit(1);
	}
	std::string str_anchor;
	std::cerr << "anchors to bin ids:\n";
	while (in >> str_anchor) {
		uint64_t anchor = str_kmer_to_uint64_t(str_anchor);
		uint64_t bin_id = MurMur64Hash{}(anchor) % n_bins;
		std::cerr << str_anchor << "\t" << bin_id << "\n";
		res[bin_id].push_back(anchor);
	}
	return res;
}

class SeparatelyOrNot {
	std::string out_path;
	bool separately;
	std::ofstream out;
public:
	SeparatelyOrNot(const std::string& out_path, bool separately) :
		out_path(out_path),
		separately(separately) {
		if (!separately) {
			out.open(out_path);
			if (!out) {
				std::cerr << "Error: cannot open file " << out_path << "\n";
				exit(1);
			}
		}
	}

	void StartBin(uint32_t bin_id) {
		if (separately) {
			std::string fname = "bin" + std::to_string(bin_id) + "." + out_path;
			out = std::ofstream(fname);
			if (!out) {
				std::cerr << "Error: cannot open file " << fname << "\n";
				exit(1);
			}
		}
	}
	std::ofstream& get_out() {
		return out;
	}
};

void process_multibin_mode(const Params& params) {
	auto bin_paths = read_bins_paths(params.input, params.n_bins);
	std::vector<std::vector<uint64_t>> anchors_bins(params.n_bins);
	auto bins_anchors = split_anchors(params.anchor_list_path, params.n_bins);
	bool accept_all = bins_anchors.empty();

	SeparatelyOrNot separately_or_not(params.output, params.separately);
	SampleNameDecoder sample_name_decoder(params.sample_names);
	for (uint32_t bin_id = 0; bin_id < params.n_bins; ++bin_id) {
		separately_or_not.StartBin(bin_id);

		if (!accept_all && bins_anchors[bin_id].empty()) {
			std::cerr << "INFO: non of specified anchors occurs in bin " << bin_id << " (" << bin_paths[bin_id]<< "). Skip reading bin content.\n";
			continue;
		}
		AcceptedAnchors accepted_anchors(accept_all ? std::vector<uint64_t>{} : bins_anchors[bin_id]);

		buffered_binary_reader in(bin_paths[bin_id]);
		if (!in) {
			std::cerr << "Error: cannot open file " << bin_paths[bin_id] << "\n";
			exit(1);
		}

		std::cerr << "Process bin " << bin_id << " (" << bin_paths[bin_id] << ")\n";

		Header header;
		header.load(in);
		header.print(std::cerr);

		Record rec;
		while (rec.load(in, header)) {
			if (accepted_anchors.IsAccepted(rec.anchor))
				rec.print(separately_or_not.get_out(), header, params.format, sample_name_decoder);
		}
	}
}

int main(int argc, char** argv) {
	if (which_bin_mode(argc, argv))
		return 0;

	auto params = read_params(argc, argv);
	params.Print(std::cerr);

	if (params.n_bins == 0)
		process_single_bin_mode(params);
	else
		process_multibin_mode(params);
}
