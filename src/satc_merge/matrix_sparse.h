#pragma once

#include "matrix_1d.h"
#include <vector>
#include <string>
#include <iostream>

namespace refresh {

template<typename I, typename T, unsigned ORDER> class matrix_sparse
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
	};

private:
	std::map<elem_t<I, ORDER>, T> data;
	size_t n_rows;
	size_t n_cols;

public:
	typedef typename std::map<elem_t<I, ORDER>, T>::iterator iterator;
	typedef typename std::map<elem_t<I, ORDER>, T>::const_iterator const_iterator;

	matrix_sparse() : n_rows(0), n_cols(0)
	{}

	matrix_sparse(size_t n_rows, size_t n_cols) : n_rows(n_rows), n_cols(n_cols)
	{}

	matrix_sparse(const matrix_sparse<I, T, ORDER>& x)
	{
		data = x.data;
		n_rows = x.n_rows;
		n_cols = x.n_cols;
	}

	matrix_sparse(matrix_sparse<I, T, ORDER>&& x)
	{
		data = move(x.data);
		n_rows = x.n_rows;
		n_cols = x.n_cols;

		x.data.clear();
		x.n_rows = 0;
		x.n_cols = 0;
	}

	matrix_sparse<I, T, ORDER>& operator=(const matrix_sparse<I, T, ORDER>& x)
	{
		if (this != &x)
		{
			data = x.data;
			n_rows = x.n_rows;
			n_cols = x.n_cols;
		}

		return *this;
	}

	matrix_sparse<I, T, ORDER>& operator=(matrix_sparse<I, T, ORDER>&& x)
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

	~matrix_sparse()
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

	T& operator()(size_t id_row, size_t id_col)
	{
		return data[elem_t<I, ORDER>(id_row, id_col)];
	}

	T operator()(size_t id_row, size_t id_col) const
	{
		return data[elem_t<I, ORDER>(id_row, id_col)];
	}

	T value(size_t id_row, size_t id_col) const
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

	const_iterator find(size_t id_row, size_t id_col)	const
	{
		return data.find(elem_t<I, ORDER>(id_row, id_col));
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

	void remove_zeros()
	{
		for (auto p = data.begin(); p != data.end();)
			if (p->second == (T)0)
				p = data.erase(p);
			else
				++p;
	}

	matrix_sparse<I, T, ORDER>& operator-=(const matrix_sparse<I, T, ORDER>& x)
	{
		for (auto& e : x.data)
			data[e.first] -= e.second;

		remove_zeros();

		return *this;
	}

	template<typename X>
	matrix_sparse<I, T, ORDER> compact(std::vector<X>& preserve_rows, std::vector<X>& preserve_cols) const 
	{
		matrix_sparse<I, T, ORDER> res(preserve_rows.size(), preserve_cols.size());

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
matrix_sparse<I, T, ORDER> operator-(const matrix_sparse<I, T, ORDER>& x, const matrix_sparse<I, T, ORDER>& y)
{
	//assert(x.rows() == y.rows && x.cols() == y.cols(), "operator- x.rows() == y.rows && x.cols() == y.cols()");

	matrix_sparse<I, T, ORDER> r(x);

	r -= y;

	return r;
}

template<typename I, typename T, unsigned ORDER>
matrix_1d<T> operator*(const matrix_1d<T>& x, const matrix_sparse<I, T, ORDER>& y)
{
	matrix_1d<T> res(y.cols(), 0);

	for (auto p = y.begin(); p != y.end(); ++p)
		res(p->first.col) += x(p->first.row) * p->second;

	return res;
}

template<typename I, typename T, unsigned ORDER>
matrix_1d<T> operator*(const matrix_sparse<I, T, ORDER>& x, const matrix_1d<T>& y)
{
	matrix_1d<T> res(x.rows(), 0);

	for (auto p = x.begin(); p != x.end(); ++p)
		res(p->first.row) += p->second * y(p->first.col);
	
	return res;
}

} // namespace refresh
