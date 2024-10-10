#include "../common/types/satc_data.h"
#include "../common/version.h"
#include "../common/accepted_anchors.h"
#include <refresh/hash_tables/lib/murmur_hash.h>
#include <string>
#include <iostream>
#include <vector>
#include <cstdint>
#include <optional>
#include <unordered_set>

#include "../common/kmc_api/kmer_defs.h"

struct Params {
	std::string input;
	std::string output;
	std::string anchor_list_path;
	std::string sample_names;
	uint32_t n_bins = 0;
	bool separately = false;
	bool binary = false;

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
		oss << "binary        : " << std::boolalpha << binary << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_dump\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <input> <output>\n";
		std::cerr << "or\n";
		std::cerr << "\t" << prog_name << " --which-bin --n_bins <int> <anchor>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input> - path to binary file in satc format\n"
			<< "    <output> - output path\n";
		std::cerr
			<< "Options:\n"
			<< "    --anchor_list <path>  - path to text file containing anchors separated by whitespaces (or tsv file with header containing column named 'anchor'), only anchors from this file will be dumped\n"
			<< "    --sample_names <path> - path for decode sample id, each line should contain <sample_name> <sample_id>\n"
			<< "    --format <string>     - output format, available options: satc, splash (default: satc)\n"
			<< "    --n_bins <int>        - if set to value different than 0 the input is interpreted as a list of bins (each bin in separate line, first line is bin_0, second line is bin_1, etc. (in case of ill-formed input results will be incorrect))\n"
			<< "    --separately          - if set with n_bins != 0 output param will be treated as suffix name and there will be output for each bin\n"
			<< "    --binary              - if set the output will be in binary satc format instead of text (satc_dump my be used as a anchor-based filter)\n";
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
	uint64_t hash = refresh::MurMur64Hash{}(anchor);
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
		else if(param == "--binary")
			res.binary = true;
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


//CRTP
template<typename impl_t, typename out_stream_t>
class writer_base {
	static_assert(std::is_same_v<out_stream_t, std::ofstream> ||
		std::is_same_v<out_stream_t, buffered_binary_writer>);
protected:
	out_stream_t out;
public:
	writer_base(const std::string& path) :
		out(path) {
		if (!out) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
	}
	void store_header(const Header& header) {
		static_cast<impl_t*>(this)->store_header(header);
	}
	void store_rec(const Record& rec, const Header& header, RecFmt format, const SampleNameDecoder& sample_name_decoder = SampleNameDecoder("")) {
		static_cast<impl_t*>(this)->store_rec(rec, header, format, sample_name_decoder);
	}
};

class writer_txt : public writer_base<writer_txt, std::ofstream>
{
public:
	writer_txt(const std::string& path) :
		writer_base<writer_txt, std::ofstream>(path) {
		
	}

	void store_header(const Header& header) {
		header.print(std::cerr);
	}

	void store_rec(Record& rec, const Header& header, RecFmt format, const SampleNameDecoder& sample_name_decoder = SampleNameDecoder("")) {
		rec.print(out, header, format, sample_name_decoder);
	}
};

class writer_binary : public writer_base<writer_binary, buffered_binary_writer>
{
public:
	writer_binary(const std::string& path) :
		writer_base<writer_binary, buffered_binary_writer>(path) {

	}

	void store_header(const Header& header) {
		header.serialize(out);
	}

	void store_rec(Record& rec, const Header& header, RecFmt format, const SampleNameDecoder& sample_name_decoder) {
		rec.serialize(out, header);
	}
};

template<typename writer_t>
void process_single_bin_mode_impl(const Params& params) {
	buffered_binary_reader in(params.input);
	if (!in) {
		std::cerr << "Error: cannot open file " << params.input << "\n";
		exit(1);
	}
	writer_t writer(params.output);
	
	Header header;
	header.load(in);

	writer.store_header(header);
	
	AcceptedAnchors accepted_anchors(params.anchor_list_path);
	Record rec;
	SampleNameDecoder sample_name_decoder(params.sample_names);
	while (rec.load(in, header)) {
		if (accepted_anchors.IsAccepted(rec.anchor))
			writer.store_rec(rec, header, params.format, sample_name_decoder);
	}
}

void process_single_bin_mode(const Params& params) {

	if(params.binary)
		process_single_bin_mode_impl<writer_binary>(params);
	else
		process_single_bin_mode_impl<writer_txt>(params);
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
	ReadAnchorsFromPlainOrDSV(in, path,
		[&res, &n_bins](uint64_t anchor)
	{
		uint64_t bin_id = refresh::MurMur64Hash{}(anchor) % n_bins;
		res[bin_id].push_back(anchor);
	});
	return res;
}


template<typename writer_t>
class SeparatelyOrNot {
	static_assert(std::is_same_v<writer_t, writer_binary> ||
		std::is_same_v<writer_t, writer_txt>);

	std::string out_path;
	bool separately;
	std::optional<writer_t> out; //no default ctor for buffered_binary_writer
	std::optional<Header> first_header;
public:
	SeparatelyOrNot(const std::string& out_path, bool separately) :
		out_path(out_path),
		separately(separately) {
		if (!separately) {
			out = writer_t(out_path);
		}
	}

	void StartBin(uint32_t bin_id) {
		if (separately) {
			std::string fname = "bin" + std::to_string(bin_id) + "." + out_path;
			out = writer_t(fname);
		}
	}
	void store_header(const Header& header) {
		if (separately)
			out->store_header(header);
		else {
			if (!first_header.has_value()) {
				first_header = header;
				out->store_header(header);
			}
			else if constexpr (std::is_same_v<writer_t, writer_binary>) {
				if (header != first_header) {
					std::cerr << "Error: input bins have different headers -> cannot be stored as a single output binary file\n";
					exit(1);
				}
			}
			else if constexpr (std::is_same_v<writer_t, writer_txt>) {
				out->store_header(header);
			}
		}
	}

	void store_rec(Record& rec, const Header& header, RecFmt format, const SampleNameDecoder& sample_name_decoder) {
		out->store_rec(rec, header, format, sample_name_decoder);
	}
};

template<typename separately_or_not_t>
void process_multibin_mode_impl(const Params& params) {

	static_assert(std::is_same_v<separately_or_not_t, SeparatelyOrNot<writer_txt>> ||
		std::is_same_v<separately_or_not_t, SeparatelyOrNot<writer_binary>>);

	auto bin_paths = read_bins_paths(params.input, params.n_bins);
	auto bins_anchors = split_anchors(params.anchor_list_path, params.n_bins);
	bool accept_all = bins_anchors.empty();

	separately_or_not_t separately_or_not(params.output, params.separately);

	SampleNameDecoder sample_name_decoder(params.sample_names);
	for (uint32_t bin_id = 0; bin_id < params.n_bins; ++bin_id) {
		separately_or_not.StartBin(bin_id);

		if (!accept_all && bins_anchors[bin_id].empty()) {
			std::cerr << "INFO: none of specified anchors occurs in bin " << bin_id << " (" << bin_paths[bin_id]<< "). Skip reading bin content.\n";
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
		separately_or_not.store_header(header);

		Record rec;
		while (rec.load(in, header)) {
			if (accepted_anchors.IsAccepted(rec.anchor))
				separately_or_not.store_rec(rec, header, params.format, sample_name_decoder);
		}
	}
}
void process_multibin_mode(const Params& params) {
	if (params.binary)
		process_multibin_mode_impl<SeparatelyOrNot<writer_binary>>(params);
	else
		process_multibin_mode_impl<SeparatelyOrNot<writer_txt>>(params);
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
