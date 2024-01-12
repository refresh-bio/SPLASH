#include "limit_mode.h"
#include <iostream>
#include <fstream>
#include "../common/version.h"
#include "../common/keep_n_largests.h"

namespace limit_mode {
	std::string to_string(Select select) {
		switch (select)
		{
		case Select::Highest:
			return "highest";
		case Select::Lowest:
			return "lowest";
		default:
			std::cerr << "Error: unhandled select " << (int)select << "\n";
			exit(1);
			break;
		}
	}

	Select select_from_string(const std::string& select) {
		if (select == "highest") {
			return Select::Highest;
		}
		else if (select == "lowest") {
			return Select::Lowest;
		}

		std::cerr << "Error: unknown select " << select << "\n";
		exit(1);
	}

	void Params::Print(std::ostream& oss) const {
		oss << "n                              : " << n << "\n";
		oss << "col_name                       : " << col_name << "\n";
		oss << "input                          : " << input_path << "\n";
		oss << "output                         : " << output_path << "\n";

		oss << "select                         : " << to_string(select) << "\n";
		oss << "separator                      : " << separator_to_string(separator) << "\n";
		oss << "sort                           : " << std::boolalpha << sort << "\n";
		oss << "sort_order                     : " << to_string(sort_order) << "\n";
		oss << "value_type                     : " << to_string(value_type) << "\n";
	}


	void Params::Usage(char* prog_name) {
		std::cerr << "dsv_manip sort mode\n";
		SPLASH_VER_PRINT(std::cerr);
		std::cerr << "Usage: \n\t" << prog_name << " limit [options] <n> <col_name> <input> <output>\n";
		std::cerr
			<< "Positional parameters:\n"
			<< "    <n>                  - number of records kept in the output\n"
			<< "    <col_name>           - name of column which values are to be sorted\n"
			<< "    <input>              - path to the input file (empty or `-` means stdin)\n"
			<< "    <output>             - path to the output file (empty or `-` means stdout)\n"
			;
		std::cerr
			<< "Options:\n"
			<< "    --select <lowest|highest>         - keep rows with lowest or highest values (default: highest)\n"
			<< "    --sep                             - separator (default: TAB)\n"
			<< "    --sort                            - sort the result\n"
			<< "    --sort_order <asc|desc>           - sort order (appliable with --sort) (default: asc)\n"
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

			if (param == "--select") {
				res.select = select_from_string(argv[++i]);
			}
			else if (param == "--sep") {
				res.separator = separator_from_string(argv[++i]);
			}
			else if (param == "--sort") {
				res.sort = true;
			}
			else if (param == "--sort_order") {
				res.sort_order = sort_order_from_string(argv[++i]);
			}
			else if (param == "--value_type") {
				res.value_type = value_type_from_string(argv[++i]);
			}
		}

		if (i >= argc) {
			std::cerr << "Error: n missing\n";
			exit(1);
		}
		res.n = std::stoull(std::string(argv[i++]));

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

	template<typename pair_t, typename PRED>
	void run_impl(const Params& params) {
		KeepNLargests<pair_t, PRED> keep_n_largests(params.n);
		using T = decltype(pair_t::first);

		std::ifstream file_in;
		std::istream& in = select_input(params.input_path, file_in);

		std::ofstream file_out;
		std::ostream& out = select_output(params.output_path, file_out);

		std::string header;
		if (!getline(in, header)) {
			std::cerr << "Error: cannot read header line\n";
			exit(1);
		}
		auto id = get_id_of_col_name(header, params.col_name, params.separator);
		std::cerr << "col id: " << id << "\n";

		std::string line;
		std::cerr << "Load the data and select...";
		while (std::getline(in, line)) {
			auto val = get_val_at<T>(line, id, params.separator);
			keep_n_largests.Add(std::make_pair(val, move(line)));
		}
		std::cerr << "\nDone.\n";
		
		std::vector<pair_t> res;
		if (params.sort) {
			std::cerr << "Sorting...";
			if (params.sort_order == SortOrder::asc)
				keep_n_largests.StealSorted(res, std::less<>{});
			else if (params.sort_order == SortOrder::desc)
				keep_n_largests.StealSorted(res, std::greater<>{});
			else {
				std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << ":" << __LINE__ << "\n";
				exit(1);
			}
			std::cerr << "\nDone.\n";
		}
		else {
			keep_n_largests.Steal(res);
		}
		
		out << header << "\n";
		for (const auto& x : res) {
			out << x.second << "\n";
		}
		
	}

	template<typename T>
	void run_impl(const Params& params) {
		using pair_t = std::pair<T, std::string>;
		
		switch (params.select)
		{
		case Select::Highest:
			run_impl<pair_t, std::greater<pair_t>>(params);
			break;
		case Select::Lowest:
			run_impl<pair_t, std::less<pair_t>>(params);
			break;
		default:
			std::cerr << "Error: this should never happen, please contact authors: " << __FILE__ << ":" << __LINE__ << "\n";
			exit(1);
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