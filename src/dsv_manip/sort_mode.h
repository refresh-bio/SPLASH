#ifndef _SORT_MODE_H
#define _SORT_MODE_H
#include <string>
#include "dsv_common.h"

namespace sort_mode {
	
	struct Params {
		char separator = '\t';
		SortOrder sort_order = SortOrder::asc;
		ValueType value_type = ValueType::_double;
		std::string col_name;
		std::string input_path; // empty or `-` for stdin
		std::string output_path; // empty or - for stdout
		
		void Print(std::ostream& oss) const;

		static void Usage(char* prog_name);
	};

	Params read_params(int argc, char** argv);

	void run(const Params& params);
}

#endif // !_SORT_MODE_H

