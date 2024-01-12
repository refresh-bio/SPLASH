#pragma once
#include <cstddef>

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

	/*	static inline uint32_t murmur_32_scramble(uint32_t k) {
			k *= 0xcc9e2d51;
			k = (k << 15) | (k >> 17);
			k *= 0x1b873593;
			return k;
		}
		uint32_t operator()(uint32_t key)
		{
			uint32_t h = 0;
			uint32_t k = key;

			key += sizeof(uint32_t);
			h ^= murmur_32_scramble(k);
			h = (h << 13) | (h >> 19);
			h = h * 5 + 0xe6546b64;

			h ^= murmur_32_scramble(k);

			h ^= 4;
			h ^= h >> 16;
			h *= 0x85ebca6b;
			h ^= h >> 13;
			h *= 0xc2b2ae35;
			h ^= h >> 16;

			return h;
		}*/
};

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

		/*		h = (~h) + (h << 21); // key = (key << 21) - key - 1;
				h = h ^ (h >> 24);
				h = (h + (h << 3)) + (h << 8); // key * 265
				h = h ^ (h >> 14);
				h = (h + (h << 2)) + (h << 4); // key * 21
				h = h ^ (h >> 28);
				h = h + (h << 31);
				return h;*/
	}
};
