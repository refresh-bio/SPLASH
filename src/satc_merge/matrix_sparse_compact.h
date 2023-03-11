#pragma once

#include "matrix_1d.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

namespace refresh {

template<typename I, typename T, unsigned ORDER> class matrix_sparse_compact
{
public:
	template<typename _I, unsigned _ORDER> struct elem_t
	{
		_I row;
		_I col;

		elem_t(_I row = 0, _I col = 0) : row(row), col(col)
		{}

		friend bool operator<(const elem_t<_I, _ORDER>& x, const elem_t<_I, _ORDER>& y)
		{
			if constexpr (_ORDER == matrix_row_major)
			{
				if (x.row != y.row)
					return x.row < y.row;
				return x.col < y.col;
			}
			else
			{
				if (x.col != y.col)
					return x.col < y.col;
				return x.row < y.row;
			}
		}

		friend bool operator==(const elem_t<_I, _ORDER>& x, const elem_t<_I, _ORDER>& y)
		{
			return x.row == y.row && x.col == y.col;
		}
	};

private:
	using data_t = std::vector<std::pair<elem_t<I, ORDER>, T>>;

	data_t data;
	size_t n_rows;
	size_t n_cols;

	size_t insert_last_into_sorted()
	{
		if (data.size() <= 1)
			return 0;

		int64_t i;
		auto elem = data.back();

		for (i = (int64_t) data.size() - 2; i >= 0; --i)
			if (elem.first < data[i].first)
				data[i + 1] = data[i];
			else
				break;

		data[i+1] = elem;

		return (size_t) (i+1);
	}

public:
	typedef typename data_t::iterator iterator;
	typedef typename data_t::const_iterator const_iterator;

	matrix_sparse_compact() : n_rows(0), n_cols(0)
	{}

	matrix_sparse_compact(size_t n_rows, size_t n_cols) : n_rows(n_rows), n_cols(n_cols)
	{}

	matrix_sparse_compact(const matrix_sparse<I, T, ORDER>& x)
	{
		data = x.data;
		n_rows = x.n_rows;
		n_cols = x.n_cols;
	}

	matrix_sparse_compact(matrix_sparse_compact<I, T, ORDER>&& x)
	{
		data = move(x.data);
		n_rows = x.n_rows;
		n_cols = x.n_cols;

		x.data.clear();
		x.n_rows = 0;
		x.n_cols = 0;
	}

	matrix_sparse_compact<I, T, ORDER>& operator=(const matrix_sparse_compact<I, T, ORDER>& x)
	{
		if (this != &x)
		{
			data = x.data;
			n_rows = x.n_rows;
			n_cols = x.n_cols;
		}

		return *this;
	}

	matrix_sparse_compact<I, T, ORDER>& operator=(matrix_sparse_compact<I, T, ORDER>&& x)
	{
		if (this != &x)
		{
			data = move(x.data);
			n_rows = x.n_rows;
			n_cols = x.n_cols;

			x.data.clear();
			x.n_rows = 0;
			x.n_cols = 0;
		}

		return *this;
	}

	~matrix_sparse_compact()
	{}

	void clear()
	{
		data.clear();
		n_rows = 0;
		n_cols = 0;
	}

	void resize(size_t _n_rows, size_t _n_cols)
	{
		data.clear();
		n_rows = _n_rows;
		n_cols = _n_cols;
	}

	void reserve(size_t _to_reserve)
	{
		data.reserve(_to_reserve);
	}

	void shrink_to_fit()
	{
		data.shrink_to_fit();
	}

	template<typename Iter>
	void assign(Iter i_begin, Iter i_end)
	{
		data.assign(i_begin, i_end);
		std::sort(data.begin(), data.end());
	}

	// Insert data to matrix without preserving ordering
	void insert_unsafe(I id_row, I id_col, const T value)
	{
		data.emplace_back(elem_t<I, ORDER>(id_row, id_col), value);
	}

	// Fix matrix ordering
	void fix()
	{
		std::sort(data.begin(), data.end());
	}

	T& operator()(I id_row, I id_col)
	{
		auto p = find(id_row, id_col);

		if (p != data.end())
			return p->second;

		data.emplace_back(elem_t<I, ORDER>(id_row, id_col), (T)0);
		auto i = insert_last_into_sorted();

		return data[i].second;
	}

/*	T operator()(I id_row, I id_col) const
	{
//		return data[elem_t<I, ORDER>(id_row, id_col)];

		return *find(id_row, id_col);
	}*/

	T value(I id_row, I id_col) const
	{
		auto p = find(id_row, id_col);

		if (p != data.end())
			return p->second;
		else
			return (T)0;
	}

	size_t rows()	const
	{
		return n_rows;
	}

	size_t cols()	const
	{
		return n_cols;
	}

	size_t data_size() const
	{
		return data.size();
	}

	const_iterator begin()	const
	{
		return data.begin();
	}

	const_iterator end()	const
	{
		return data.end();
	}

	T sum()	const
	{
		T r = (T)0;

		for (auto& x : data)
			r += x.second;

		return r;
	}

	bool exists(size_t id_row, size_t id_col)	const
	{
		return find(id_row, id_col) != data.end();
	}

	const_iterator find(I id_row, I id_col)	const
	{
		elem_t<I, ORDER> elem{ id_row, id_col };

		if (data.empty())
			return data.end();

		// When matrix is filled often the accessed element could be the last one7
		if (data.back().first == elem)
			return --data.end();

		std::pair<elem_t<I, ORDER>, T> item{ elem, (T)0 };

		auto p = std::lower_bound(data.begin(), data.end(), item, [](const auto& x, const auto& y) {
			return x.first < y.first;
			});

		if (p == data.end())
			return p;

		if (p->first == elem)
			return p;

		return data.end();
	}

	iterator find(I id_row, I id_col)
	{
		if (data.empty())
			return data.end();

		elem_t<I, ORDER> elem{ id_row, id_col };

		// When matrix is filled often the accessed element could be the last one7
		if (data.back().first == elem)
			return --data.end();
		if (data.back().first < elem)
			return data.end();

		std::pair<elem_t<I, ORDER>, T> item{ elem, (T)0 };

		auto p = std::lower_bound(data.begin(), data.end(), item, [](const auto& x, const auto &y) {
			return x.first < y.first;
			});

		if (p == data.end())
			return p;

		if (p->first == elem)
			return p;

		return data.end();
	}

	void set_all_to_zero()
	{
		data.clear();
	}

	matrix_1d<T> get_row_sums()	const
	{
		matrix_1d<T> vec(n_rows, 0);

		for (const auto& x : data)
			vec(x.first.row) += x.second;

		return vec;
	}

	matrix_1d<T> get_col_sums()	const
	{
		matrix_1d<T> vec(n_cols, 0);

		for (const auto& x : data)
			vec(x.first.col) += x.second;

		return vec;
	}

	void get_row_col_sums(matrix_1d<T>& vec_rows, matrix_1d<T>& vec_cols)	const
	{
		vec_rows.clear();
		vec_cols.clear();

		vec_rows.resize(n_rows, 0);
		vec_cols.resize(n_cols, 0);

		for (const auto& x : data)
		{
			vec_rows(x.first.row) += x.second;
			vec_cols(x.first.col) += x.second;
		}
	}

	matrix_sparse_compact<I, T, ORDER>& operator-=(matrix_sparse_compact<I, T, ORDER>& x)
	{
		data_t new_elems;

		auto p_src1 = data.begin();
		auto p_src2 = x.data.begin();
		auto p_dest = data.begin();

		auto p_end1 = data.end();
		auto p_end2 = x.data.end();

		while (p_src1 != p_end1 && p_src2 != p_end2)
		{
			auto elem = *p_src1;

			if (p_src1->first == p_src2->first)
			{
				elem.second -= p_src2->second;
				if(elem.second != (T)0)
					*p_dest++ = elem;
				++p_src1;
				++p_src2;
			}
			else if (p_src1->first < p_src2->first)
			{
				*p_dest++ = elem;
				++p_src1;
			}
			else
			{
				new_elems.emplace_back(*p_src2);
				++p_src2;
			}
		}

		if (p_src1 != p_end1)
			while (p_src1 != p_end1)
				*p_dest++ = *p_src1++;
		else
			new_elems.insert(new_elems.end(), p_src2, p_end2);

		data.erase(p_dest, p_end1);

		if (!new_elems.empty())
		{
			for (auto& elem : new_elems)
				elem.second = -elem.second;

			data.insert(data.end(), new_elems.begin(), new_elems.end());
			
			// TODO: consider using inplace_merge here
			std::sort(data.begin(), data.end());
		}

		return *this;
	}

	template<typename X>
	matrix_sparse_compact<I, T, ORDER> compact(std::vector<X>& preserve_rows, std::vector<X>& preserve_cols) const
	{
		// TODO: Consider in-place relabeling
		matrix_sparse_compact<I, T, ORDER> res(preserve_rows.size(), preserve_cols.size());

		res.reserve(data.size());

		auto rows_begin = preserve_rows.begin();
		auto rows_end = preserve_rows.end();
		auto cols_begin = preserve_cols.begin();
		auto cols_end = preserve_cols.end();

		for (auto& x : data)
		{
			auto p_row = std::lower_bound(rows_begin, rows_end, x.first.row);
			if (p_row == rows_end || *p_row != x.first.row)
				continue;

			auto p_col = std::lower_bound(cols_begin, cols_end, x.first.col);
			if (p_col == cols_end || *p_col != x.first.col)
				continue;

			I new_row = p_row - rows_begin;
			I new_col = p_col - cols_begin;

			res(new_row, new_col) = x.second;
		}

		res.shrink_to_fit();

		return res;
	}

	template<typename X>
	matrix_sparse_compact<I, T, ORDER> compact(std::vector<X>& preserve_rows) const
	{
		// TODO: Consider in-place relabeling
		matrix_sparse_compact<I, T, ORDER> res(preserve_rows.size(), cols());

		res.reserve(data.size());

		auto rows_begin = preserve_rows.begin();
		auto rows_end = preserve_rows.end();

		for (auto& x : data)
		{
			auto p_row = std::lower_bound(rows_begin, rows_end, x.first.row);
			if (p_row == rows_end || *p_row != x.first.row)
				continue;

			I new_row = p_row - rows_begin;

			res(new_row, x.first.col) = x.second;
		}

		res.shrink_to_fit();

		return res;
	}

	void print(const std::string& name) const
	{
#ifdef REFRESH_MATRIX_ALLOW_PRINT
		std::cout << name << ":\n[";
		for (size_t i = 0; i < n_rows; ++i)
		{
			std::cout << "\n[ ";
			for (size_t j = 0; j < n_cols; ++j)
			{
				std::cout << value(i, j);
				std::cout << "  ";
			}
			std::cout << "]";
		}
		std::cout << "]\n";
#endif
	}

};

template<typename I, typename T, unsigned ORDER>
matrix_sparse_compact<I, T, ORDER> operator-(const matrix_sparse_compact<I, T, ORDER>& x, const matrix_sparse_compact<I, T, ORDER>& y)
{
	matrix_sparse_compact<I, T, ORDER> r(x);

	r -= y;

	return r;
}

template<typename I, typename T, unsigned ORDER>
matrix_1d<T> operator*(const matrix_1d<T>& x, const matrix_sparse_compact<I, T, ORDER>& y)
{
	matrix_1d<T> res(y.cols(), 0);

	for (auto p = y.begin(); p != y.end(); ++p)
		res(p->first.col) += x(p->first.row) * p->second;

	return res;
}

template<typename I, typename T, unsigned ORDER>
matrix_1d<T> operator*(const matrix_sparse_compact<I, T, ORDER>& x, const matrix_1d<T>& y)
{
	matrix_1d<T> res(x.rows(), 0);

	for (auto p = x.begin(); p != x.end(); ++p)
		res(p->first.row) += p->second * y(p->first.col);
	
	return res;
}

} // namespace refresh
