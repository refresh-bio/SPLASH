#ifndef _BLOOM_SET_H
#define _BLOOM_SET_H

#include <algorithm>
#include <utility>
#include <cmath>

namespace refresh
{
// **********************************************************************************
template<typename Key_t,
	typename Hash_t = std::hash<Key_t>,
	unsigned NO_HASHES = 3>
class bloom_set {

public:
	using key_type = Key_t;
	using hasher = Hash_t;

private:
	Hash_t hash;

	uint64_t* arr = nullptr;
	uint64_t* raw_arr = nullptr;

	size_t no_elements;
	size_t allocated;
	size_t mask;
	uint32_t mask_shift;

	static size_t normalize_size(size_t size, double expected_fp_rate)
	{
		size_t expected_no_items = size;

		size *= NO_HASHES;
		size *= 2;

		if (size & (size - 1))
		{
			while (size & (size - 1))
				size &= size - 1;
			size *= 2;
		}

		size = std::max((size_t)256, size);

		while (est_fp_rate(expected_no_items, size) > expected_fp_rate && size <= (1ull << 60))
			size *= 2;

		return size;
	}

	static double est_fp_rate(size_t no_items, size_t size_in_bits)
	{
		return std::pow(1 - std::exp(-(double)NO_HASHES * no_items / (double)size_in_bits), NO_HASHES);
	}

	void do_allocation()
	{
		raw_arr = new uint64_t[allocated / 64 + 7];
		arr = raw_arr;
		while (((uint64_t)arr) % 64 != 0)
			++arr;
	}

	void allocate(size_t size, double expected_fp_rate)
	{
		if (raw_arr)
			delete[] raw_arr;
		no_elements = 0;

		allocated = normalize_size(size, expected_fp_rate);

		do_allocation();

		std::fill_n(arr, allocated / 64, 0ull);

		mask_shift = 6 * NO_HASHES;
		mask = (allocated / 64 - 1) << mask_shift;
	}

	void insert_impl(const size_t h)
	{
	//	uint64_t h = mmh(x);
		size_t pos = (h & mask) >> mask_shift;
		uint64_t or_mask = 0;

		if constexpr (NO_HASHES == 2)
			or_mask = (1ull << (h & 63)) | (1ull << ((h >> 6) & 63));
		else if constexpr (NO_HASHES == 3)
			or_mask = (1ull << (h & 63)) | (1ull << ((h >> 6) & 63)) | (1ull << ((h >> 12) & 63));

		arr[pos] |= or_mask;

		++no_elements;
	}

	bool check_impl(const size_t h) const
	{
		size_t pos = (h & mask) >> mask_shift;

		if constexpr (NO_HASHES == 2)
			return (arr[pos] & (1ull << (h & 63))) && (arr[pos] & (1ull << ((h >> 6) & 63)));
		else if constexpr (NO_HASHES == 3)
			return (arr[pos] & (1ull << (h & 63))) && (arr[pos] & (1ull << ((h >> 6) & 63))) && (arr[pos] & (1ull << ((h >> 12) & 63)));
	}

public:
	bloom_set(size_t size = 64, double expected_fp_rate = 0.05)
	{
		static_assert(NO_HASHES == 2 || NO_HASHES == 3, "No. of hashes could be 2 or 3");
		allocate(size, std::max(expected_fp_rate, 0.001));
	}

	bloom_set(const bloom_set<Key_t, Hash_t, NO_HASHES>& rhs)
	{
		hash = rhs.hash;

		no_elements = rhs.no_elements;
		allocated = rhs.allocated;
		mask = rhs.mask;
		mask_shift = rhs.mask_shift;

		do_allocation();

		std::copy_n(rhs.arr, allocated, arr);
	}

	bloom_set(bloom_set<Key_t, Hash_t, NO_HASHES>&& rhs) noexcept
	{
		hash = rhs.hash;

		no_elements = rhs.no_elements;
		allocated = rhs.allocated;
		mask = rhs.mask;
		mask_shift = rhs.mask_shift;

		arr = rhs.arr;
		raw_arr = rhs.raw_arr;

		rhs.arr = nullptr;
		rhs.raw_arr = nullptr;
	}

	bloom_set& operator=(bloom_set<Key_t, Hash_t, NO_HASHES>& rhs)
	{
		if (this != &rhs)
		{
			if (raw_arr)
				delete[] raw_arr;

			hash = rhs.hash;

			no_elements = rhs.no_elements;
			allocated = rhs.allocated;
			mask = rhs.mask;
			mask_shift = rhs.mask_shift;

			do_allocation();

			std::copy_n(rhs.arr, allocated, arr);
		}

		return *this;
	}

	bloom_set& operator=(bloom_set<Key_t, Hash_t, NO_HASHES>&& rhs) noexcept
	{
		if (this != &rhs)
		{
			if (raw_arr)
				delete[] raw_arr;

			hash = rhs.hash;

			no_elements = rhs.no_elements;
			allocated = rhs.allocated;
			mask = rhs.mask;
			mask_shift = rhs.mask_shift;

			arr = rhs.arr;
			raw_arr = rhs.raw_arr;

			rhs.arr = nullptr;
			rhs.raw_arr = nullptr;
		}

		return *this;
	}

	~bloom_set()
	{
		if (raw_arr)
			delete[] raw_arr;
	}

	// 0 - no change
	void resize(size_t size, double expected_fp_rate = 0.05)
	{
		allocate(size, std::max(expected_fp_rate, 0.001));
	}

	template<typename Iter>
	void insert(Iter begin, Iter end)
	{
		for (auto p = begin; p != end; ++p)
			insert_impl(hash(*p));
	}

	void insert(const Key_t &x)
	{
		insert_impl(hash(x));
	}

	bool check(const Key_t &x) const
	{
		return check_impl(hash(x));
	}

	bool check(const Key_t &x, size_t h) const
	{
		return check_impl(h);
	}

	double filling_factor()
	{
		return (double)NO_HASHES * no_elements / allocated;
	}

	static size_t capacity(double expected_fp_rate = 0.05)
	{
		int bits_for_sub_hashes = 6 * NO_HASHES;

		for (size_t i = 1; i < 64; ++i)
		{
			size_t req_size_in_bits = normalize_size(1ull << i, expected_fp_rate);
			size_t req_no_words = req_size_in_bits / 64;

			if (req_no_words > (1ull << (64 - bits_for_sub_hashes)))
				return 1ull << (i-1);
		}

		return 0;
	}

	size_t mem_usage()
	{
		return allocated / 8;
	}

	double fp_rate()
	{
		return est_fp_rate(no_elements, allocated);
	}
};
}

#endif