#include "csv.h"
#include <algorithm>
#include <numeric>

namespace refresh
{
// ************************************************************************************
// csv_basic
// ************************************************************************************
// ************************************************************************************
void csv_basic::split(const std::string& str, std::vector<std::string>& parts)
{
	parts.clear();

	std::string s;

	for (auto c : str)
	{
		if (in_sep_map[c])
		{
			parts.emplace_back(s);
			s.clear();
		}
		else
			s.push_back(c);
	}

	if (!s.empty())
		parts.emplace_back(s);
}


// ************************************************************************************
// csv_file
// ************************************************************************************

// ************************************************************************************
std::string csv_file::join(size_t row_id)
{
	std::string str;

	if (row_id == std::numeric_limits<size_t>::max())		// Header
		for (size_t i = 0; i < data.size(); ++i)
			str += data[i].first + out_sep;
	else
		for (size_t i = 0; i < data.size(); ++i)
			str += data[i].second[row_id] + out_sep;

	str.pop_back();

	return str;
}
 
// ************************************************************************************
bool csv_file::load(const std::string& fn)
{
	std::ifstream ifs(fn);

	if (ifs.bad())
	{
//		cerr << "Cannot open file: " << fn << endl;
		return false;
	}

	std::vector<std::string> parts;
	std::string line;

	getline(ifs, line);
	split(line, parts);

	data.clear();
	data.reserve(parts.size());

	for (auto& str : parts)
		data.emplace_back(str, std::vector<std::string>());

	size_t no_parts = data.size();

	while (!ifs.eof())
	{
		getline(ifs, line);
		if (ifs.eof())
			break;

		split(line, parts);
		if (parts.size() != no_parts)
			return false;

		for (size_t i = 0; i < no_parts; ++i)
			data[i].second.emplace_back(parts[i]);
	}

	return true;
}

// ************************************************************************************
void csv_file::initialize(const csv_file& csv)
{
	data.clear();
	data.shrink_to_fit();

	data.reserve(csv.data.size());

	for (const auto& x : csv.data)
		data.emplace_back(x.first, std::vector<std::string>());
}

// ************************************************************************************
void csv_file::mark(const size_t row_id)
{
	if (marked.size() != no_rows())
	{
		marked.clear();
		marked.resize(no_rows(), false);
	}

	marked[row_id] = true;
}

// ************************************************************************************
void csv_file::remove_marked()
{
	size_t i_orig, i_after;
	size_t nc = no_cols();
	size_t nr = no_rows();

	for (i_orig = i_after = 0; i_orig < nr; ++i_orig)
	{
		if (marked[i_orig])
			continue;
		if (i_orig != i_after)
			for (size_t j = 0; j < nc; ++j)
				data[j].second[i_after] = data[j].second[i_orig];
		++i_after;
	}

	for (size_t j = 0; j < nc; ++j)
	{
		data[j].second.resize(i_after, "");
		data[j].second.shrink_to_fit();
	}

	marked.clear();
	marked.shrink_to_fit();
}

// ************************************************************************************
size_t csv_file::no_cols() const
{
	return data.size();
}

// ************************************************************************************
size_t csv_file::no_rows() const
{
	if (data.empty())
		return 0;
	return data.front().second.size();
}

// ************************************************************************************
bool csv_file::save(const std::string& fn)
{
	std::ofstream ofs(fn);

	if (ofs.bad())
		return false;

	ofs << join() << std::endl;

	size_t nr = no_rows();

	for (size_t i = 0; i < nr; ++i)
		ofs << join(i) << std::endl;

	return true;
}

// ************************************************************************************
std::vector<std::string> csv_file::header() const
{
	std::vector<std::string> names;

	names.reserve(data.size());

	for (const auto& x : data)
		names.emplace_back(x.first);

	return names;
}

// ************************************************************************************
void csv_file::rename_col(size_t col_id, const std::string& new_name)
{
	data[col_id].first = new_name;
}

// ************************************************************************************
bool csv_file::merge(const csv_file& csv)
{
	if (data.empty())
	{
		data = csv.data;
		return true;
	}

	auto h1 = header();
	auto h2 = csv.header();

	if (h1 != h2)
		return false;

	for (size_t i = 0; i < data.size(); ++i)
		data[i].second.insert(data[i].second.end(), csv.data[i].second.begin(), csv.data[i].second.end());

	return true;
}

// ************************************************************************************
bool csv_file::merge(csv_file&& csv)
{
	if (data.empty())
	{
		data = move(csv.data);
		csv.clear();

		return true;
	}

	auto h1 = header();
	auto h2 = csv.header();

	if (h1 != h2)
		return false;

	size_t new_no_rows = no_rows() + csv.no_rows();
	size_t csv_no_rows = csv.no_rows();

	for (size_t i = 0; i < data.size(); ++i)
	{
		auto& col = data[i].second;
		auto& csv_col = csv.data[i].second;

		col.reserve(new_no_rows);

		for (size_t j = 0; j < csv_no_rows; ++j)
			col.emplace_back(move(csv_col[j]));

		csv_col.clear();
		csv_col.shrink_to_fit();
	}

	csv.clear();

	return true;
}

// ************************************************************************************
void csv_file::clear()
{
	data.clear();
	data.shrink_to_fit();
}

// ************************************************************************************
int csv_file::col_id(const std::string& col_name)
{
	for (int i = 0; i < (int)data.size(); ++i)
		if (data[i].first == col_name)
			return i;

	return -1;
}

// ************************************************************************************
std::vector<std::string>& csv_file::col(size_t col_id)
{
	return data[col_id].second;
}

// ************************************************************************************
std::vector<std::string>& csv_file::col(const std::string& col_name)
{
	return col(col_id(col_name));
}

// ************************************************************************************
void csv_file::remove_col(const size_t col_id)
{
	data.erase(data.begin() + col_id);
}

// ************************************************************************************
void csv_file::remove_col(const std::string& col_name)
{
	int c_id = col_id(col_name);

	if(c_id >= 0)
		remove_col(c_id);
}

// ************************************************************************************
std::vector<std::string> csv_file::row_view(const size_t row_id, const std::initializer_list<size_t> col_ids) const
{
	std::vector<std::string> row;

	if (col_ids.size() == 0)
	{
		row.reserve(data.size());

		for (size_t i = 0; i < data.size(); ++i)
			row.emplace_back(data[i].second[row_id]);
	}
	else
	{
		row.reserve(col_ids.size());

		for (auto id : col_ids)
			row.emplace_back(data[id].second[row_id]);
	}

	return row;
}

// ************************************************************************************
bool csv_file::append_row(std::vector<std::string> parts)
{
	if (parts.size() != data.size())
		return false;

	for (size_t i = 0; i < parts.size(); ++i)
		data[i].second.emplace_back(parts[i]);

	return true;
}

// ************************************************************************************
void csv_file::sort()
{
	auto nr = no_rows();

	if (nr < 2)
		return;

	std::vector<size_t> order(nr);

	std::iota(order.begin(), order.end(), 0);

	std::stable_sort(order.begin(), order.end(), [&](const auto x, const auto y)
		{
			for (size_t i = 0; i < data.size(); ++i)
				if (data[i].second[x] != data[i].second[y])
					return data[i].second[x] < data[i].second[y];

			return x < y;
		});

	std::vector<std::string> tmp_col;

	for (size_t i = 0; i < data.size(); ++i)
	{
		auto& col = data[i].second;

		tmp_col = move(col);
		col.clear();
		col.reserve(nr);

		for (size_t j = 0; j < nr; ++j)
			col.emplace_back(move(tmp_col[order[j]]));
	}

	if (!marked.empty())
	{
		std::vector<bool> tmp_marked;

		tmp_marked.reserve(nr);

		for (size_t i = 0; i < nr; ++i)
			tmp_marked.emplace_back(marked[order[i]]);

		marked = move(tmp_marked);
	}
}

// ************************************************************************************
// csv_stream
// ************************************************************************************

// ************************************************************************************
// csv_istream
// ************************************************************************************
bool csv_istream::open(const std::string& filename)
{
	ifs.open(filename, std::ios_base::binary);

	ifs.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);

	load_header();

	return ifs.good();
}

// ************************************************************************************
bool csv_istream::read_line()
{
	cur_line.clear();

	std::getline(ifs, buffer);

	return !buffer.empty();
}

// ************************************************************************************
bool csv_istream::load_header() 
{
	header.clear();

	getline(ifs, buffer);
	split(buffer, header);

	return !header.empty();
}

// ************************************************************************************
bool csv_istream::get_record(std::vector<std::string>& record)
{
	if (!ifs.good())
		return false;

	getline(ifs, buffer);
	split(buffer, record);

	return record.size() == header.size();
}

// ************************************************************************************
// csv_ostream
// ************************************************************************************
bool csv_ostream::open(const std::string& filename)
{
	ofs.open(filename, std::ios_base::binary);

	ofs.rdbuf()->pubsetbuf(io_buffer, io_buffer_size);

	header_stored = false;

	return ofs.good();
}

// ************************************************************************************
bool csv_ostream::set_header(const std::vector<std::string>& col_names)
{
	if (!ofs.is_open())
		return false;

	header = col_names;

	return !header.empty();
}

// ************************************************************************************
bool csv_ostream::add_record(const std::vector<std::string>& record)
{
	if (!ofs.good())
		return false;

	if (!header_stored)
	{
		join(header);
		ofs.write(buffer.c_str(), buffer.size());
		header_stored = true;
	}

	if (header.size() != record.size())
		return false;

	join(record);

	ofs.write(buffer.c_str(), buffer.size());

	return ofs.good();
}

// ************************************************************************************
void csv_ostream::join(const std::vector<std::string>& cur_line)
{
	buffer.clear();

	for (const auto& x : cur_line)
	{
		buffer += x;
		buffer.push_back(out_sep);
	}

	buffer.pop_back();
	buffer.push_back('\n');
}

}