#include "../common/satc_data.h"
#include "../common/version.h"
#include "../common/accepted_anchors.h"
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

	RecFmt format = RecFmt::SATC;
	void Print(std::ostream& oss) const
	{
		oss << "input         : " << input << "\n";
		oss << "output        : " << output << "\n";
		oss << "anchor_list   : " << anchor_list_path << "\n";
		oss << "sample_names  : " << sample_names << "\n";
		oss << "format        : " << RecFmtConv::to_string(format) << "\n";
	}

	static void Usage(char* prog_name) {
		std::cerr << "satc_dump\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " [options] <input> <output>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input> - path to binary file in satc format\n"
			<< "    <output> - output path\n";
		std::cerr 
			<< "Options:\n"
			<< "    --anchor_list <path>  - path to text file containing anchors separated by whitespaces, only anchors from this file will be dumped\n"
			<< "    --sample_names <path> - path for decode sample id, each line should contain <sample_name> <sample_id>\n"
			<< "    --format <string>     - output format, available options: satc, splash (default: satc)\n";
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

		if (param == "--anchor_list")
			res.anchor_list_path = argv[++i];
		else if(param == "--sample_names")
			res.sample_names = argv[++i];
		else if (param == "--format")
			res.format = RecFmtConv::from_string(argv[++i]);
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

int main(int argc, char** argv) {

	auto params = read_params(argc, argv);
	params.Print(std::cerr);

	buffered_binary_reader in(params.input);
	if (!in) {
		std::cerr << "Error: cannot open file " << params.input << "\n";
		return 1;
	}

	std::ofstream out(params.output);
	if (!out) {
		std::cerr << "Error: cannot open file " << params.output << "\n";
		return 1;
	}

	Header header;
	header.load(in);
	header.print(std::cerr);

	AcceptedAnchors accepted_anchors(params.anchor_list_path);
	Record rec;
	SampleNameDecoder sample_name_decoder(params.sample_names);
	while (rec.load(in, header)) {
		if(accepted_anchors.IsAccepted(rec.anchor))
			rec.print(out, header, params.format, sample_name_decoder);
	}

}
