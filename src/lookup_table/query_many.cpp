#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "query_many.h"
#include "../common/version.h"

namespace query_many_mode {

	void Params::Print(std::ostream& oss) const {
		query_cfg.Print(oss);
	}

	void Params::Usage(char* prog_name) {
		std::cerr << "lookup_table (query lookup table with many files)\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " query_many <input> <lookup>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <input>              - path to the file where each line defines sequence files to query, each line follows the format of regular query mode\n"
			<< "    <lookup>             - path to the index build using `build` mode\n"
			;

		std::cerr
			<< "Options:\n"
			<< "    --n_threads <int> - number of threads\n"
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

		for (; i < argc; ++i) {
			if (argv[i][0] != '-')
				break;

			std::string param = argv[i];
			if (param == "--n_threads") {
				std::string tmp = argv[++i];
				res.query_cfg.n_threads = std::stoull(tmp);
			}
		}

		if (i >= argc) {
			std::cerr << "Error: input missing\n";
			exit(1);
		}
		res.input = argv[i++];
		
		if (i >= argc) {
			std::cerr << "Error: lookup missing\n";
			exit(1);
		}
		res.query_cfg.lookup = argv[i++];

		std::ifstream in(res.input);
		if (!in) {
			std::cerr << "Error: cannot open file" << res.input << "\n";
			exit(1);
		}

		std::string line;
		while (std::getline(in, line)) {
			std::istringstream iss(line);
			std::vector<std::string> config { std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{ } } ;
			res.query_cfg.input_cfg.emplace_back(config);
		}

		return res;
	}

	void run(const Params& params) {
		query_common::run(params.query_cfg);
	}
}