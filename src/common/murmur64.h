#pragma once
#include <cstddef>

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
