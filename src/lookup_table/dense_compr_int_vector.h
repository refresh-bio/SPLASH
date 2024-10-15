#ifndef _DENSE_COMPRESSED_INT_VECTOR_H
#define _DENSE_COMPRESSED_INT_VECTOR_H
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>

//value_type is a type that allows to store values from iv
template<typename value_type>
class dense_compr_int_vector
{
	//simple class that construct some data structure to fast query for existence
	//not sure what is the best, for my use case I will start with unordered_set
	template<typename Key_t, typename Val_t>
	class check_dense
	{
		std::unordered_map<Key_t, Val_t> s;
	public:

		void add(const Key_t& key, const Val_t& val)
		{
			s.emplace(key, val);
		}

		bool contains(const Key_t& key, Val_t& val)
		{
			auto it = s.find(key);
			if (it != s.end())
			{
				val = it->second;
				return true;
			}
			return false;
		}
	};

	//using rank_type = sdsl::rank_support_v5<>;
	using rank_type = sdsl::rank_support_v<>;

	std::vector<value_type> dense_values;
	sdsl::int_vector<> sparse_values;
	sdsl::int_vector<> indexes_to_dense_values;
	rank_type rs_bv;

	sdsl::bit_vector bv;

	size_t tot_size_bytes{};

	std::vector<std::pair<size_t, value_type>> get_freq_and_values(const sdsl::int_vector<>& iv)
	{
		//TODO: potential consideration is that the range of values in iv may be quite large making this histogram large, so may be better to use some dictionary data structure in such a case
		//std::cerr << "Computing histogram...";
		std::vector<size_t> histogram;
		for (auto val : iv)
		{
			if (val >= histogram.size())
			{
				histogram.reserve(val * 2 + 1);
				histogram.resize(val + 1);
			}
			++histogram[val];
		}
		//std::cerr << "Done.\n";

		//std::cerr << "Converting histogram..\n";
		//                    freq    value
		std::vector<std::pair<size_t, value_type>> freq_and_values;
		
		size_t tot_sum{};

		for (size_t val = 0; val < histogram.size(); ++val)
		{
			if (histogram[val])
			{
				freq_and_values.emplace_back(histogram[val], val);

				tot_sum += histogram[val];
			}
		}

		//double best_possible_compression = 0;
		//
		//for (const auto& x : freq_and_values)
		//	best_possible_compression += -log2((double)x.first / tot_sum) * x.first;
		//
		//std::cerr << "Best possible compression (bits?): " << (size_t)best_possible_compression << "\n";
		//std::cerr << "Best possible compression (bytes?): " << (size_t)(best_possible_compression / 8) << "\n";

		//std::cerr << "Done.\n";

		//std::cerr << "Sorting...\n";
		std::sort(freq_and_values.begin(), freq_and_values.end(), std::greater<>{});

		//std::cerr << "Done.\n";

		return freq_and_values;
	}

	void accumulate_frequencies(std::vector<std::pair<size_t, value_type>>& freq_and_values)
	{
		//std::cerr << "Accumulate frequencies...\n";
		size_t prev{};
		for (auto& x : freq_and_values)
			prev = x.first += prev;

		//std::cerr << "Done.\n";
	}

	size_t find_best_b(const std::vector<std::pair<size_t, value_type>>& freq_and_values, size_t n, size_t p)
	{
		//n - number of elements
		//p - bits/elem
		//b - bits for dense element - to be found
		//x - number of dense elements

		//the total cost (in bits) of compressed variant is:
		// n         - bit vector to distinguish dense and sparse values in the original int_vector
		// rs        - cost of rank support
		// (x-n) * p - cost to store sparse values
		// x * p     - cost to store dense values (actually it works a little faster
		//             when std::vector is used instead of int_vector of width b), so this
		//             cost is about x * sizeof(value_type)
		// x * 2^b   - cost of indexes of dense values

				//in the loop below we find b that minimizes total cost described above
		//I skip cost of rank_support because its not very important
		//but in general could be included
		//and may be considered
		size_t best_size_bits = std::numeric_limits<size_t>::max();
		size_t best_b{};

		//could start with b = 0, but this would requier special handling
		//in this case we have only one dense value, so the index to them is not needed
		for (size_t b = 1; b < p; ++b)
		{
			size_t pow_2_p = 1ull << b;

			if (pow_2_p - 1 >= freq_and_values.size())
				break;
			size_t x = freq_and_values[pow_2_p - 1].first;

			size_t bv_size_bits = n;

			size_t dense_values_size_bits = (1ull << b) * sizeof(value_type);

			size_t sparse_values_size_bits = (n - x) * p;

			size_t indexes_to_dense_values_size_bits = x * b;

			size_t total_size_bits = bv_size_bits + dense_values_size_bits + sparse_values_size_bits + indexes_to_dense_values_size_bits;

			//std::cerr << b << "\t" << total_size_bits << "\n";

			if (total_size_bits < best_size_bits)
			{
				best_size_bits = total_size_bits;
				best_b = b;
			}
		}
		//at this point we could verify if best_size_bits is better than plain representation
		//and make a decision

		assert(best_size_bits != std::numeric_limits<size_t>::max());

		return best_b;
	}

public:
	//not exact, but should be pretty close
	size_t get_tot_size_bytes() const
	{
		return tot_size_bytes;
	}

	dense_compr_int_vector() = default;
	dense_compr_int_vector(const sdsl::int_vector<>& iv)
	{
		static_assert(std::is_unsigned_v<value_type>);
		if (8 * sizeof(value_type) < iv.width())
		{
			std::cerr << "Error: the width ( " << iv.width() << ") of a given int_vector is larger than may be stored using " << 8 * sizeof(value_type) << " bits\n";
			exit(1);
		}

		//no data
		if (iv.empty())
			return;

		auto freq_and_values = get_freq_and_values(iv);

		
		accumulate_frequencies(freq_and_values);
		
		size_t n = iv.size();
		size_t p = iv.width();

		size_t b = find_best_b(freq_and_values, n, p);

		//std::cerr << "best b: " << b << "\n";

		size_t pow_2_p = 1ull << b;

		//number of dense elements
		size_t x = freq_and_values[pow_2_p - 1].first;

		//std::cerr << "Collect dense values...\n";

		dense_values = std::vector<value_type>(1ull << b);

		check_dense<value_type, uint32_t> cd;

		for (size_t i = 0; i < pow_2_p; ++i) {
			auto key = freq_and_values[i].second;
			cd.add(key, i);
			dense_values[i] = key;
		}

		//std::cerr << "Done.\n";

		//std::cerr << "Build...\n";
		sparse_values = sdsl::int_vector<>(n - x, 0, p);
		indexes_to_dense_values = sdsl::int_vector<>(x, 0, b);
		size_t pos_indexes_to_dense_values{};
		size_t pos_sparse_values{};

		bv = sdsl::bit_vector(n);

		for (size_t i = 0; i < n; ++i)
		{
			value_type val = iv[i];
			uint32_t idx;
			if (cd.contains(val, idx))
			{
				bv[i] = true;
				indexes_to_dense_values[pos_indexes_to_dense_values++] = idx;
			}
			else
			{
				sparse_values[pos_sparse_values++] = val;
			}
		}

		//std::cerr << "Done.\n";

		//std::cerr << "Build rank support for bit_vector\n";

		rs_bv = rank_type(&bv);

		//std::cerr << "Done\n";

		//std::cerr << "bv size in bytes: " << sdsl::size_in_bytes(bv) << "\n";
		//std::cerr << "sparse_values size in bytes: " << sdsl::size_in_bytes(sparse_values) << "\n";
		//std::cerr << "indexes_to_dense_values bytes: " << sdsl::size_in_bytes(indexes_to_dense_values) << "\n";
		//std::cerr << "dense values size in bytes: " << sdsl::size_in_bytes(dense_values) << "\n";
		//std::cerr << "rs_bv bytes: " << sdsl::size_in_bytes(rs_bv) << "\n";
		tot_size_bytes = sdsl::size_in_bytes(bv) + sdsl::size_in_bytes(dense_values) + 
			sdsl::size_in_bytes(sparse_values) + sdsl::size_in_bytes(indexes_to_dense_values) + 
			sdsl::size_in_bytes(rs_bv);

		//std::cerr << "tot size in bytes: " << tot_size_bytes << "\n";
	}

	value_type get(size_t idx)
	{
		if (bv[idx])
		{
			auto rank = rs_bv.rank(idx);
			return dense_values[indexes_to_dense_values[rank]];
		}
		else
		{
			auto rank_0 = idx - rs_bv.rank(idx);
			return sparse_values[rank_0];
		}
	}

	void serialize(std::ostream& out) const
	{
		bv.serialize(out);
		indexes_to_dense_values.serialize(out);

		//create sdsl int_vector for dense values (this should be small)
		//because this "sdsl::serialize(dense_values, out)" didn't work I dont know why,
		//but dont have time to investigate 
		sdsl::int_vector<sizeof(value_type) * 8> tmp(dense_values.size());
		std::copy(dense_values.begin(), dense_values.end(), tmp.begin());
		tmp.serialize(out);
		
		sparse_values.serialize(out);
	}

	void load(std::istream& in)
	{
		bv.load(in);
		rs_bv = rank_type(&bv);

		indexes_to_dense_values.load(in);

		sdsl::int_vector<sizeof(value_type) * 8> tmp;
		tmp.load(in);
		dense_values = std::vector<value_type>(tmp.begin(), tmp.end());

		sparse_values.load(in);
	}

	value_type operator[](size_t idx)
	{
		return get(idx);
	}
};



#endif // _DENSE_COMPRESSED_INT_VECTOR_H
