#pragma once
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include "matrix.h"

namespace detail {
	//my iteartor to not generate the range...
	struct MyIt {

		using iterator_category = std::random_access_iterator_tag;
		using value_type = size_t;
		using reference = value_type&;
		using difference_type = size_t;

		size_t v;
		MyIt(size_t v) : v(v) {}


		reference operator*() {
			return v;
		}
		MyIt& operator++() {
			++v;
			return *this;
		}
		difference_type operator-(const MyIt& rhs) const {
			return v - rhs.v;
		}

		bool operator==(const MyIt& rhs) {
			return v == rhs.v;
		}
	};
}

template<> struct std::iterator_traits<detail::MyIt>
{
	using iterator_category = typename detail::MyIt::iterator_category;
	using value_type = typename detail::MyIt::value_type;
	using difference_type = typename detail::MyIt::difference_type;
};

template<typename OutIter, typename RNG>
void my_sample(size_t start, size_t end, OutIter dest, size_t N, RNG&& eng) {
	std::sample(detail::MyIt(start), detail::MyIt(end), dest, N, eng);
}

template<typename RNG>
//refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major> get_train_mtx_2(const refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major>& X, double train_fraction, RNG&& eng) 
refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> get_train_mtx_2(const refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major>& X, double train_fraction, RNG&& eng) 
{
	//refresh::matrix_sparse<uint32_t, double, refresh::matrix_col_major> res(X.rows(), X.cols());
	refresh::matrix_sparse_compact<uint32_t, double, refresh::matrix_col_major> res(X.rows(), X.cols());

	res.reserve(X.data_size());

	std::vector<double> cumsum(X.rows());

	auto p = X.begin();

	for (size_t col_id = 0; col_id < X.cols(); ++col_id) {
		fill(cumsum.begin(), cumsum.end(), 0.0);

		for (; p != X.end() && p->first.col == col_id; ++p)
			cumsum[p->first.row] = p->second;
		std::partial_sum(cumsum.begin(), cumsum.end(), cumsum.begin());
		
		size_t to_select = cumsum.back() * train_fraction;
		std::vector<size_t> indices(to_select);
		my_sample(0, cumsum.back(), indices.begin(), to_select, eng);

		for (size_t i = 0, j = 0; i < cumsum.size() && j < indices.size();)
			if (indices[j] < cumsum[i]) {
				++res(i, col_id);
				++j;
			}
			else
				++i;
	}

	res.shrink_to_fit();

	return res;
}