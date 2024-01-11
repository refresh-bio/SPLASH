#pragma once

#include "matrix_1d.h"
#include <string>
#include <iostream>

namespace refresh {

template<typename T, unsigned ORDER> class matrix_dense
{
	T* ptr;
	size_t n_rows;
	size_t n_cols;

	void allocate()
	{
		ptr = new T[n_rows * n_cols];
	}

	void release()
	{
		if (ptr)
		{
			delete[] ptr;
			ptr = nullptr;

			n_rows = 0;
			n_cols = 0;
		}
	}

public:
	matrix_dense() : ptr(nullptr), n_rows(0), n_cols(0) {}

	matrix_dense(size_t n_rows, size_t n_cols) : ptr(nullptr), n_rows(n_rows), n_cols(n_cols)
	{
		allocate();
	}

	matrix_dense(const matrix_dense<T, ORDER>& x)
	{
		n_rows = x.n_rows;
		n_cols = x.n_cols;

		allocate();

		std::memcpy(ptr, x.ptr, n_rows * n_cols * sizeof(T));
	}

	matrix_dense(matrix_dense<T, ORDER>&& x)
	{
		n_rows = x.n_rows;
		n_cols = x.n_cols;
		ptr = x.ptr;

		x.n_rows = 0;
		x.n_cols = 0;
		x.ptr = nullptr;
	}

	~matrix_dense()
	{
		release();
	}

	matrix_dense<T, ORDER>& operator=(const matrix_dense<T, ORDER>& x)
	{
		if (this != &x)
		{
			release();
			n_rows = x.n_rows;
			n_cols = x.n_cols;
			allocate();

			std::memcpy(ptr, x.ptr, n_rows * n_cols * sizeof(T));
		}

		return *this;
	}

	matrix_dense<T, ORDER>& operator=(matrix_dense<T, ORDER>&& x)
	{
		if (this != &x)
		{
			release();

			n_rows = x.n_rows;
			n_cols = x.n_cols;
			ptr = x.ptr;

			x.n_rows = 0;
			x.n_cols = 0;
			x.ptr = nullptr;
		}

		return *this;
	}

	void clear()
	{
		release();
	}

	size_t rows()
	{
		return n_rows;
	}

	size_t cols()
	{
		return n_cols;
	}

	size_t size()
	{
		return n_rows * n_cols;
	}

	T* data()
	{
		return ptr;
	}

	T& operator()(size_t id_row, size_t id_col)
	{
		if constexpr (ORDER == matrix_row_major)
			return ptr[id_row * n_rows + id_col];
		else
			return ptr[id_col * n_rows + id_row];
	}

	T operator()(size_t id_row, size_t id_col) const
	{
		if constexpr (ORDER == matrix_row_major)
			return ptr[id_row * n_rows + id_col];
		else
			return ptr[id_col * n_rows + id_row];
	}

	T sum()
	{
		return std::accumulate(ptr, ptr + n_rows * n_cols, (T)0);
	}

	void set_to_zero()
	{
		std::fill_n(ptr, n_rows * n_cols, (T)0);
	}

	matrix_1d<T> get_col(size_t col_id)
	{
		matrix_1d<T> vec(n_rows);

		if constexpr (ORDER == matrix_row_major)
			for (size_t i = 0; i < n_rows; ++i)
				vec(i) = ptr[col_id + i * n_cols];
		else
			for (size_t i = 0; i < n_rows; ++i)
				vec(i) = ptr[col_id * n_rows + i];

		return vec;
	}

	matrix_1d<T> get_row(size_t row_id)
	{
		matrix_1d<T> vec(n_cols);

		if constexpr (ORDER == matrix_row_major)
			for (size_t i = 0; i < n_cols; ++i)
				vec(i) = ptr[row_id * n_cols + i];
		else
			for (size_t i = 0; i < n_cols; ++i)
				vec(i) = ptr[row_id + i * n_rows];

		return vec;
	}

	matrix_1d<T> get_row_sums()
	{
		matrix_1d<T> vec(n_rows, 0);

		if constexpr (ORDER == matrix_row_major)
		{
			for (size_t i = 0; i < n_rows; ++i)
			{
				T sum = (T)0;
				auto p = ptr + i * n_cols;

				for (size_t j = 0; j < n_cols; ++j)
					sum += p[j];

				vec(i) = sum;
			}
		}
		else
		{
			auto pv = vec.data();
			auto pm = ptr;

			for (int i = 0; i < n_cols; ++i)
				for (int j = 0; j < n_rows; ++j)
					pv[j] += *pm++;
		}

		return vec;
	}

	matrix_1d<T> get_col_sums()
	{
		matrix_1d<T> vec(n_cols, 0);

		if constexpr (ORDER == matrix_row_major)
		{
			auto pv = vec.data();
			auto pm = ptr;

			for (int i = 0; i < n_rows; ++i)
				for (int j = 0; j < n_cols; ++j)
					pv[j] += *pm++;
		}
		else
			for (size_t i = 0; i < n_cols; ++i)
			{
				T sum = (T) 0;
				auto p = ptr + i * n_rows;

				for (size_t j = 0; j < n_rows; ++j)
					sum += p[j];

				vec(i) = sum;
			}

		return vec;
	}
};

} // namespace refresh
