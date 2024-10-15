#ifndef _MURMUR_HASH_H
#define _MURMUR_HASH_H

#include <cstdint>
#include <cstddef>

namespace refresh {
	// **********************************************************************************
	struct MurMur32Hash
	{
		std::size_t operator()(uint32_t h) const noexcept
		{
			h ^= h >> 16;
			h *= 0x85ebca6b;
			h ^= h >> 13;
			h *= 0xc2b2ae35;
			h ^= h >> 16;

			return h;
		}
	};

	// **********************************************************************************
	// MurMurHash3
	struct MurMur64Hash
	{
		std::size_t operator()(size_t h) const noexcept
		{
			h ^= h >> 33;
			h *= 0xff51afd7ed558ccdL;
			h ^= h >> 33;
			h *= 0xc4ceb9fe1a85ec53L;
			h ^= h >> 33;

			return h;
		}
	};

	// **********************************************************************************
	// MurMurHash3 for pair
	struct MurMurPair64Hash
	{
		std::size_t operator()(const std::pair<std::size_t, std::size_t>& x) const noexcept
		{
			std::size_t h = x.first;

			h ^= h >> 33;
			h *= 0xff51afd7ed558ccdL;
			h ^= h >> 33;
			h *= 0xc4ceb9fe1a85ec53L;
			h ^= h >> 33;

			h ^= x.second;

			h ^= h >> 33;
			h *= 0xff51afd7ed558ccdL;
			h ^= h >> 33;
			h *= 0xc4ceb9fe1a85ec53L;
			h ^= h >> 33;

			return h;
		}
	};

	// **********************************************************************************
	/// MurMurHash3 for strings (simple implementation)
	struct MurMurStringsHash
	{
		// Based on https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
	private:
		static uint64_t load64(const char*& p)
		{
			uint64_t x = (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);
			x <<= 8;		x += (uint64_t)(*p++);

			return x;
		}

		static inline uint64_t rotl64(uint64_t x, int8_t r)
		{
			return (x << r) | (x >> (64 - r));
		}

		static inline uint64_t fmix64(uint64_t k)
		{
			k ^= k >> 33;
			k *= 0xff51afd7ed558ccdull;
			k ^= k >> 33;
			k *= 0xc4ceb9fe1a85ec53ull;
			k ^= k >> 33;

			return k;
		}

	public:
		std::size_t operator()(const std::string& s) const
		{
			uint64_t h1 = 0;
			uint64_t h2 = 0;

			//		const std::size_t n_blocks = s.size() / 16;

			const uint64_t c1 = 0x87c37b91114253d5ull;
			const uint64_t c2 = 0x4cf5ad432745937full;

			const char* data = s.c_str();

			for (std::size_t i = 0; i < s.size() / 16; i++)
			{
				uint64_t k1 = load64(data);
				uint64_t k2 = load64(data);

				k1 *= c1; k1 = rotl64(k1, 31); k1 *= c2; h1 ^= k1;

				h1 = rotl64(h1, 27); h1 += h2; h1 = h1 * 5 + 0x52dce729;

				k2 *= c2; k2 = rotl64(k2, 33); k2 *= c1; h2 ^= k2;

				h2 = rotl64(h2, 31); h2 += h1; h2 = h2 * 5 + 0x38495ab5;
			}

			std::size_t tail = s.size() % 16;

			uint64_t k1 = 0;
			uint64_t k2 = 0;

			switch (tail & 15)
			{
			case 15: k2 ^= ((uint64_t)data[14]) << 48; [[fallthrough]];
			case 14: k2 ^= ((uint64_t)data[13]) << 40; [[fallthrough]];
			case 13: k2 ^= ((uint64_t)data[12]) << 32; [[fallthrough]];
			case 12: k2 ^= ((uint64_t)data[11]) << 24; [[fallthrough]];
			case 11: k2 ^= ((uint64_t)data[10]) << 16; [[fallthrough]];
			case 10: k2 ^= ((uint64_t)data[9]) << 8; [[fallthrough]];
			case  9: k2 ^= ((uint64_t)data[8]) << 0;
				k2 *= c2; k2 = rotl64(k2, 33); k2 *= c1; h2 ^= k2;
				[[fallthrough]];
			case  8: k1 ^= ((uint64_t)data[7]) << 56; [[fallthrough]];
			case  7: k1 ^= ((uint64_t)data[6]) << 48; [[fallthrough]];
			case  6: k1 ^= ((uint64_t)data[5]) << 40; [[fallthrough]];
			case  5: k1 ^= ((uint64_t)data[4]) << 32; [[fallthrough]];
			case  4: k1 ^= ((uint64_t)data[3]) << 24; [[fallthrough]];
			case  3: k1 ^= ((uint64_t)data[2]) << 16; [[fallthrough]];
			case  2: k1 ^= ((uint64_t)data[1]) << 8; [[fallthrough]];
			case  1: k1 ^= ((uint64_t)data[0]) << 0;
				k1 *= c1; k1 = rotl64(k1, 31); k1 *= c2; h1 ^= k1;
			};

			h1 ^= (uint64_t)s.size(); h2 ^= (uint64_t)s.size();

			h1 += h2;
			h2 += h1;

			h1 = fmix64(h1);
			h2 = fmix64(h2);

			h1 += h2;
			h2 += h1;

			return h1 ^ h2;
		}
	};
}
#endif