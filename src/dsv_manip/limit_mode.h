#ifndef _LIMIT_MODE_H
#define _LIMIT_MODE_H
#include "dsv_common.h"
#include <string>

namespace limit_mode {
	enum class Select { Lowest, Highest };
	struct Params {

		Select select = Select::Highest;
		char separator = '\t';
		bool sort = false;
		SortOrder sort_order = SortOrder::asc;
		ValueType value_type = ValueType::_double;
		
		size_t n;
		std::string col_name;
		std::string input_path; // empty or `-` for stdin
		std::string output_path; // empty or - for stdout
		
		void Print(std::ostream& oss) const;

		static void Usage(char* prog_name);
	};

	Params read_params(int argc, char** argv);

	void run(const Params& params);
}

#endif // !_LIMIT_MODE_H

