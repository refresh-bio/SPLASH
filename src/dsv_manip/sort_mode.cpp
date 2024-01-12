#include "sort_mode.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "../common/version.h"

namespace sort_mode {

	void Params::Print(std::ostream& oss) const {
		oss << "col_name                       : " << col_name << "\n";
		oss << "input                          : " << input_path << "\n";
		oss << "output                         : " << output_path << "\n";

		oss << "separator                      : " << separator_to_string(separator) << "\n";
		oss << "sort_order                     : " << to_string(sort_order) << "\n";
		oss << "value_type                     : " << to_string(value_type) << "\n";
	}

	void Params::Usage(char* prog_name) {
		std::cerr << "dsv_manip sort mode\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " sort [options] <col_name> <input> <output>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <col_name>           - name of column which values are to be sorted\n"
			<< "    <input>              - path to the input file (empty or `-` means stdin)\n"
			<< "    <output>             - path to the output file (empty or `-` means stdout)\n"
			;
		std::cerr
			<< "Options:\n"
			<< "    --sep                             - separator (default: TAB)\n"
			<< "    --sort_order <asc|desc>           - sort order (default: asc)\n"
			<< "    --value_type <int|double|string>  - the type of values stored in column with col_name (default: double)\n"
			;
	}

	Params read_params(int argc, char** argv) {
		Params res;
		if (argc == 2) {
			Params::Usage(argv[0]);
			exit(0);
		}
		int i = 2;
		for (; i < argc; ++i)
		{
			if (argv[i][0] != '-')
				break;

			std::string param = argv[i];

			if (param == "--sep") {
				res.separator = separator_from_string(argv[++i]);
			}
			else if (param == "--sort_order") {
				res.sort_order = sort_order_from_string(argv[++i]);
			}
			else if (param == "--value_type") {
				res.value_type = value_type_from_string(argv[++i]);
			}
		}

		if (i >= argc) {
			std::cerr << "Error: col_name missing\n";
			exit(1);
		}
		res.col_name = argv[i++];

		if (i < argc) {
			res.input_path = argv[i++];
		}
		else {
			res.input_path = "-";
		}

		if (i < argc) {
			res.output_path = argv[i++];
		}
		else {
			res.output_path = "-";
		}

		return res;
	}

	

	template<typename T>
	void run_impl(const Params& params) {
		std::ifstream file_in;
		std::istream& in = select_input(params.input_path, file_in);

		std::ofstream file_out;
		std::ostream& out = select_output(params.output_path, file_out);

		std::vector<std::pair<T, std::string>> data;

		std::string header;
		if (!getline(in, header)) {
			std::cerr << "Error: cannot read header line\n";
			exit(1);
		}
		auto id = get_id_of_col_name(header, params.col_name, params.separator);
		std::cerr << "col id: " << id << "\n";

		std::string line;
		std::cerr << "Loading the data...";
		while (std::getline(in, line)) {
			auto val = get_val_at<T>(line, id, params.separator);
			data.emplace_back(val, move(line));
		}

		std::cerr << "\nDone.\n";
		std::cerr << "Sorting...";
		if (params.sort_order == SortOrder::asc)
			std::sort(data.begin(), data.end(), std::less<>{});
		else if (params.sort_order == SortOrder::desc)
			std::sort(data.begin(), data.end(), std::greater<>{});
		else {
			std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << ":" << __LINE__ << "\n";
			exit(1);
		}
		std::cerr << "\nDone.\n";

		out << header << "\n";
		for (const auto& x : data) {
			out << x.second << "\n";
		}
	}

	void run(const Params& params) {
		switch (params.value_type)
		{
		case ValueType::_double:
			run_impl<double>(params);
			break;
		case ValueType::_int:
			run_impl<int64_t>(params);
			break;
		case ValueType::_string:
			run_impl<std::string>(params);
			break;
		default:
			std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << ":" << __LINE__ << "\n";
			exit(1);
		}
	}
}