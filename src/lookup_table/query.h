#ifndef QUERY_H
#define QUERY_H
#include <string>
#include <cinttypes>
#include <vector>
#include "query_common.h"

namespace query_mode {
	using namespace query_common;

	struct Params {
		QueryCfg query_cfg;

		void Print(std::ostream& oss) const;

		static void Usage(char* prog_name);
	};

	Params read_params(int argc, char** argv);

	void run(const Params& params);
}

#endif // !QUERY_H
