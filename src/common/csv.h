#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <limits>
#include <cinttypes>
#include <algorithm>

namespace refresh
{
	// ************************************************************************************
	// ************************************************************************************
	class csv_basic
	{
	protected:
		char out_sep;
		std::string in_seps;
		bool in_sep_map[256];

		void split(const std::string& str, std::vector<std::string>& parts);

		// ************************************************************************************
		std::string my_to_string(double x)
		{
			char buf[32];
			sprintf(buf, "%.12g", x);
			return std::string(buf);
		}

		// ************************************************************************************
		std::string my_to_string(std::string x)
		{
			return x;
		}

		// ************************************************************************************
		double my_to_double(const std::string& str)
		{
			double x;
			if (sscanf(str.c_str(), "%lf", &x) != 1)
			{
				std::cerr << "sscanf error\n";
				exit(1);
			}
			return x;
		}

		// ************************************************************************************
		std::string my_to_string(int64_t x)
		{
			return std::to_string(x);
		}

		// ************************************************************************************
		void prepare_in_sep_map()
		{
			std::fill_n(in_sep_map, 256, false);
			for (unsigned char c : in_seps)
				in_sep_map[c] = true;
		}

	public:
		csv_basic(const std::string in_seps, const char out_sep) : in_seps(in_seps), out_sep(out_sep)
		{
			prepare_in_sep_map();
		}

		csv_basic(const char in_sep, const char out_sep) : out_sep(out_sep)
		{
			in_seps.clear();
			in_seps.push_back(in_sep);

			prepare_in_sep_map();
		}
	};

	// ************************************************************************************
	// ************************************************************************************
	class csv_file : public csv_basic
	{
		std::vector<std::pair<std::string, std::vector<std::string>>> data;

		std::vector<bool> marked;

		std::string join(size_t row_id = std::numeric_limits<size_t>::max());

	public:
		csv_file(const std::string in_seps = ",", const char out_sep = ',') : csv_basic(in_seps, out_sep)
		{}

		csv_file(const char in_sep = ',', const char out_sep = ',') : csv_basic(in_sep, out_sep)
		{}

		void change_separator(const std::string _in_seps, const char _out_sep)
		{
			in_seps = _in_seps;
			out_sep = _out_sep;

			prepare_in_sep_map();
		}

		size_t no_cols() const;
		size_t no_rows() const;

		bool load(const std::string& fn);
		bool save(const std::string& fn);

		void initialize(const csv_file& csv);

		std::vector<std::string> header() const;

		bool empty() const
		{
			return no_rows() == 0;
		}

		void clear();

		bool merge(const csv_file& csv);
		bool merge(csv_file&& csv);

		int col_id(const std::string& col_name);

		std::vector<std::string>& col(size_t col_id);
		std::vector<std::string>& col(const std::string& col_name);

		void rename_col(size_t col_id, const std::string& new_name);

		std::vector<std::string> copy_col_str(size_t col_id)
		{
			return data[col_id].second;
		}

		std::vector<std::string> copy_col_str(const std::string& col_name)
		{
			return data[col_id(col_name)].second;
		}

		std::vector<int64_t> copy_col_int(size_t col_id)
		{
			auto& col = data[col_id].second;
			std::vector<int64_t> res;

			res.reserve(no_rows());

			for (auto& x : col)
				res.emplace_back(stoll(x));

			return res;
		}

		std::vector<int64_t> copy_col_int(const std::string& col_name)
		{
			return copy_col_int(col_id(col_name));
		}

		std::vector<double> copy_col_double(size_t col_id)
		{
			auto& col = data[col_id].second;
			std::vector<double> res;

			res.reserve(no_rows());

			for (auto& x : col)
				res.emplace_back(my_to_double(x));

			return res;
		}

		std::vector<double> copy_col_double(const std::string& col_name)
		{
			return copy_col_double(col_id(col_name));
		}

		/*		bool replace_col(const size_t col_id, std::vector<T>& col)
				{
					if (col.size() != no_rows || col_id < no_cols())
						return false;

					auto csv_col = data[col_id].second;
					auto nr = col.size();

					for (size_t i = 0; i < nr; ++i)
						csv_col[i] = my_to_string(col[i]);

					return true;
				}

				template<typename T> bool replace_col(const std::string &col_name, std::vector<T>& col)
				{
					return replace_col(col_id(col_name));
				}*/

		template<typename T> bool insert_col(const std::string& col_name, std::vector<T>& col)
		{
			if (col.size() != no_rows() && no_rows() != 0)
				return false;

			if (col_id(col_name) >= 0)
				return false;

			data.emplace_back(col_name, std::vector<std::string>());

			auto& csv_col = data.back().second;
			auto nr = col.size();

			csv_col.reserve(nr);

			for (size_t i = 0; i < nr; ++i)
				csv_col.emplace_back(my_to_string(col[i]));

			return true;
		}

		void mark(const size_t row_id);
		void remove_marked();

		void remove_col(const size_t col_id);
		void remove_col(const std::string& col_name);

		std::vector<std::string> row_view(const size_t row_id, const std::initializer_list<size_t> col_ids) const;
		bool append_row(std::vector<std::string> parts);

		template<typename Iter>
csv_file filter(Iter begin, Iter end)
{
	csv_file csv_out(in_seps, out_sep);

	for (auto p = begin; p != end; ++p)
		csv_out.insert_col(data[*p].first, data[*p].second);

	return csv_out;
}

void sort();
	};

	// ************************************************************************************
	// ************************************************************************************
	class csv_stream : public csv_basic
	{
	protected:
		std::vector<std::string> header;
		std::string buffer;

		const size_t io_buffer_size = 16 << 20;
		char* io_buffer;

	public:
		csv_stream(const std::string in_seps, const char out_sep) : csv_basic(in_seps, out_sep)
		{
			io_buffer = new char[io_buffer_size];
		}

		csv_stream(const char in_sep, const char out_sep) : csv_basic(in_sep, out_sep)
		{
			io_buffer = new char[io_buffer_size];
		}

		~csv_stream()
		{
			delete[] io_buffer;
		}

		csv_stream(csv_stream& x) = delete;

		int col_id(const std::string& col_name)
		{
			auto p = std::find(header.begin(), header.end(), col_name);
			if (p == header.end())
				return -1;
			return p - header.begin();
		}
	};

	// ************************************************************************************
	// ************************************************************************************
	class csv_istream : public csv_stream
	{
		std::ifstream ifs;

		std::vector<std::string> cur_line;

		bool load_header();

		bool read_line();

	public:
		csv_istream(const std::string& filename = "", const std::string in_seps = ",") : csv_stream(in_seps, 0)
		{
			if (!filename.empty())
				open(filename);
		}

		csv_istream(const std::string& filename, const char in_sep = ',') : csv_stream(in_sep, 0)
		{
			if (!filename.empty())
				open(filename);
		}

		csv_istream(csv_istream& x) = delete;

		bool good()
		{
			return ifs.good();
		}

		bool is_open()
		{
			return ifs.is_open();
		}

		bool eof()
		{
			return ifs.eof();
		}

		bool open(const std::string& filename);

		std::vector<std::string>& get_header()
		{
			return header;
		}

		bool get_record(std::vector<std::string>& record);

		double to_double(const std::string& s)
		{
			return my_to_double(s);
		}
	};

	// ************************************************************************************
	// ************************************************************************************
	class csv_ostream : public csv_stream
	{
		std::ofstream ofs;
		bool header_stored;

		void join(const std::vector<std::string> &cur_line);

	public:
		csv_ostream(const std::string& filename, const char out_sep = ',') : csv_stream("", out_sep), header_stored(false)
		{
			if (!filename.empty())
				open(filename);
		}

		csv_ostream(csv_ostream& x) = delete;

		bool open(const std::string& filename);

		bool good()
		{
			return ofs.good();
		}

		bool is_open()
		{
			return ofs.is_open();
		}

		bool set_header(const std::vector<std::string>& col_names);
		bool add_record(const std::vector<std::string>& record);

		std::string to_string(const double d)
		{
			return my_to_string(d);
		}
	};
}

