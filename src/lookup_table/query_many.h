#ifndef QUERY_MANY_H
#define QUERY_MANY_H
#include "query_common.h"

namespace query_many_mode {
	using namespace query_common;

	struct Params {
		QueryCfg query_cfg;

		std::string input;

		void Print(std::ostream& oss) const;

		static void Usage(char* prog_name);
	};

	Params read_params(int argc, char** argv);

	void run(const Params& params);
}

#endif // !QUERY_MANY_H