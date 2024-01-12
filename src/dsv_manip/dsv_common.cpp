#include "dsv_common.h"
#include <string>
#include <iostream>
#include <fstream>

std::string to_string(SortOrder sort_order) {
	switch (sort_order)
	{
	case SortOrder::asc:
		return "asc";
	case SortOrder::desc:
		return "desc";
	default:
		std::cerr << "Error: unhandled sort_order " << (int)sort_order << "\n";
		exit(1);
		break;
	}
}

SortOrder sort_order_from_string(const std::string& sort_order) {
	if (sort_order == "asc") {
		return SortOrder::asc;
	}
	else if (sort_order == "desc") {
		return SortOrder::desc;
	}

	std::cerr << "Error: unknown sort_order " << sort_order << "\n";
	exit(1);
}


std::string to_string(ValueType value_type) {
	switch (value_type)
	{
	case ValueType::_double:
		return TypeName<double>::name;
	case ValueType::_int:
		return TypeName<int64_t>::name;
	case ValueType::_string:
		return TypeName<std::string>::name;
	default:
		std::cerr << "Error: unhandled value_type " << (int)value_type << "\n";
		exit(1);
		break;
	}
}

ValueType value_type_from_string(const std::string& value_type) {
	if (value_type == TypeName<double>::name) {
		return ValueType::_double;
	}
	else if (value_type == TypeName<int64_t>::name) {
		return ValueType::_int;
	}
	else if (value_type == TypeName<std::string>::name) {
		return ValueType::_string;
	}

	std::cerr << "Error: unknown value_type " << value_type << "\n";
	exit(1);
}


std::string separator_to_string(char sep) {
	if (sep == '\t')
		return "TAB";
	return std::string(1, sep);
}

char separator_from_string(const std::string& sep) {
	if (sep == "TAB")
		return '\t';
	if (sep.length() != 1) {
		std::cerr << sep << " cannot be used as separator\n";
		exit(1);
	}
	return sep[0];
}

std::istream& select_input(const std::string& path, std::ifstream& file_in) {
	if (path != "-") {
		file_in.open(path);
		if (!file_in) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		return file_in;
	}
	return std::cin;
}

std::ostream& select_output(const std::string& path, std::ofstream& file_out) {
	if (path != "-") {
		file_out.open(path);
		if (!file_out) {
			std::cerr << "Error: cannot open file " << path << "\n";
			exit(1);
		}
		return file_out;
	}
	return std::cout;
}

int get_id_of_col_name(const std::string& header, const std::string& col_name, char sep) {
	int i = 0;
	std::istringstream iss(header);
	std::string tmp;
	while (std::getline(iss, tmp, sep)) {
		if (tmp == col_name)
			return i;
		++i;
	}
	std::cerr << "Error: cannot find column " << col_name << "\n";
	exit(1);
}
