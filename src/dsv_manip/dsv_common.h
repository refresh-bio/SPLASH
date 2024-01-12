#ifndef _DSV_COMMON_H
#define _DSV_COMMON_H
#include <string>
#include <iostream>
#include <sstream>

enum class SortOrder { asc, desc };
enum class ValueType { _double, _int, _string };

template<typename T>
struct TypeName;

template<>
struct TypeName<double>
{
	inline static const std::string name = "double";
};
template<>
struct TypeName<int64_t>
{
	inline static const std::string name = "int";
};
template<>
struct TypeName<std::string>
{
	inline static const std::string name = "string";
};

std::string to_string(SortOrder sort_order);

SortOrder sort_order_from_string(const std::string& sort_order);

std::string to_string(ValueType value_type);

ValueType value_type_from_string(const std::string& value_type);

std::string separator_to_string(char sep);

char separator_from_string(const std::string& sep);

std::istream& select_input(const std::string& path, std::ifstream& file_in);

std::ostream& select_output(const std::string& path, std::ofstream& file_out);

int get_id_of_col_name(const std::string& header, const std::string& col_name, char sep);
template<typename T>
inline bool convert_to(const std::string& x, T& res) {
	std::istringstream iss(x);
	if (!(iss >> res))
		return false;
	return true;
}

template<>
inline bool convert_to<std::string>(const std::string& x, std::string& res) {
	res = x;
	return true;
}

template<typename T>
T get_val_at(const std::string& line, int col_id, char sep) {
	int i = 0;
	std::istringstream iss(line);
	std::string tmp;
	while (std::getline(iss, tmp, sep)) {
		if (i == col_id) {
			T res;
			if (!convert_to(tmp, res)) {
				std::cerr << "Error: cannot convert " << tmp << " to " << TypeName<T>::name << "\n";
				exit(1);
			}
			return res;
		}
		++i;
	}
	std::cerr << "Error: cannot find column id " << col_id << "\n";
	exit(1);
}
#endif // !_DSV_COMMON_H

