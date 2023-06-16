#pragma once
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include "matrix.h"

#include <iostream>
using std::cerr;
using std::endl;

template<typename RNG>
void my_sample(size_t max_val, std::vector<size_t>& dest, size_t to_select, RNG&& eng) {
	dest.resize(max_val);
	std::iota(dest.begin(), dest.end(), 0ull);

	for (size_t i = 0; i < to_select; ++i)
		std::swap(dest[i], dest[i + eng() % (dest.size() - i)]);

	dest.resize(to_select);
	sort(dest.begin(), dest.end());
}

template<typename RNG>
refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> get_train_mtx_2(
	const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X,
	double opt_train_fraction,
	RNG&& eng)
{
	refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> res(X.rows(), X.cols());

	res.reserve(X.data_size());

	std::vector<size_t> indices;

	auto p = X.begin();

	for (size_t col_id = 0; col_id < X.cols(); ++col_id) {	
		double tot_cs = 0;
		auto p_col_end = p;

		for (; p_col_end != X.end() && p_col_end->first.col == col_id; ++p_col_end)
			tot_cs += p_col_end->second;

		size_t to_select = tot_cs * opt_train_fraction;
		my_sample(tot_cs, indices, to_select, eng);
		indices.push_back(std::numeric_limits<size_t>::max()); // guard

		double cs = 0;

		for (int j = 0; p != p_col_end; ++p)
		{
			cs += p->second;

			for (; indices[j] < cs; ++j)
				++res(p->first.row, col_id);
		}
	}

	res.shrink_to_fit();

	return res;
}