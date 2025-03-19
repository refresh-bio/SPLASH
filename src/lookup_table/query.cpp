#include "query.h"
#include <iostream>
#include <fstream>
#include "../common/version.h"


namespace query_mode {

	void Params::Print(std::ostream& oss) const {
		query_cfg.Print(oss);
	}

	void Params::Usage(char* prog_name) {
		std::cerr << "lookup_table (query lookup table)\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " query [options] <input> <lookup> <output>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input>              - path to the file with sequences to be queried\n"
			<< "    <lookup>             - path to the index build using `build` mode\n"
			<< "    <output>             - output file, each line contains results for each sequence\n"
			;
		std::cerr
			<< "Options:\n"
			<< "    --input_fmt <string>   - input format, one of: fasta, extendors, compactors (default: fasta)\n"
			<< "    --report_fmt <string>  - format of the detailed report, one of: plain (verbose format), concise (default), ids (just ids), empty (do not report, useful when stats enabled with stats_fmt) (default: concise)\n"
			<< "    --stats_fmt <string>   - for each query report stats (#k-mers per category) one of: empty (don't print stats) or with_stats (print stats) (default: empty)\n"
			<< "    --output_fmt <string>  - output format, one of: txt (one line for each query), extendors (only for extendors input, adds additional columns to the extendors file), compactors (default: txt)\n"
			<< "    --kmer_skip <int>      - for each query sequence the next queried k-mer start on position <kmer_skip> + 1\n"
			<< "    --n_threads <int>      - number of threads\n"
			;
	}

	Params read_params(int argc, char** argv)
	{
		Params res;
		if (argc == 2) {
			Params::Usage(argv[0]);
			exit(0);
		}
		int i = 2;

		std::vector<std::string> input_cfg;

		for (; i < argc; ++i)
		{
			if (argv[i][0] != '-')
				break;

			std::string param = argv[i];

			if (param == "--input_fmt") {
				input_cfg.push_back("--input_fmt");
				input_cfg.push_back(argv[++i]);
			}
			if (param == "--report_fmt") {
				input_cfg.push_back("--report_fmt");
				input_cfg.push_back(argv[++i]);
			}
			if (param == "--stats_fmt") {
				input_cfg.push_back("--stats_fmt");
				input_cfg.push_back(argv[++i]);
			}
			if (param == "--output_fmt") {
				input_cfg.push_back("--output_fmt");
				input_cfg.push_back(argv[++i]);
			}
			if (param == "--kmer_skip") {
				input_cfg.push_back("--kmer_skip");
				input_cfg.push_back(argv[++i]);
			}
			if (param == "--n_threads") {
				res.query_cfg.n_threads = true;
			}
		}

		if (i >= argc) {
			std::cerr << "Error: input missing\n";
			exit(1);
		}
		input_cfg.push_back(argv[i++]);

		if (i >= argc) {
			std::cerr << "Error: lookup missing\n";
			exit(1);
		}
		res.query_cfg.lookup = argv[i++];

		if (i >= argc) {
			std::cerr << "Error: output missing\n";
			exit(1);
		}
		input_cfg.push_back(argv[i++]);

		res.query_cfg.input_cfg.emplace_back(input_cfg);

		return res;
	}

	void run(const Params& params) {
		query_common::run(params.query_cfg);
	}
}
