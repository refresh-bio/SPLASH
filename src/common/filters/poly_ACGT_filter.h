#ifndef _POLY_ACGT_FILTER_H
#define _POLY_ACGT_FILTER_H
#include <cinttypes>

class PolyACGTFilter
{
	uint32_t len;
	uint64_t mask;

	uint64_t patternC;
	uint64_t patternG;
	uint64_t patternT;

	uint64_t build_pattern(uint32_t symb)
	{
		uint64_t res = symb;
		for (uint32_t i = 1; i < len; ++i)
		{
			res <<= 2;
			res += symb;
		}
		return res;
	}

	bool is_poly(uint64_t part) const {
		return part == 0        || //poly(A)
			   part == patternC ||
			   part == patternG ||
			   part == patternT;
	}

public:
	PolyACGTFilter(uint32_t len) //len==0 means no filtering
		:len(len),
		mask(((1ull)<<(2*len))-1),
		patternC(build_pattern(1)),
		patternG(build_pattern(2)),
		patternT(build_pattern(3))
	{

	}
	PolyACGTFilter() :
		PolyACGTFilter(0)
	{

	}

	auto GetLen() const { return len; }
	bool IsPolyACGT(uint64_t mer, uint32_t len) const
	{
		if (this->len == 0)
			return false;
		if (this->len > len)
			return false;
		auto no_parts = len - this->len + 1;
		for (uint32_t i = 0; i < no_parts - 1; ++i)
		{
			if (is_poly(mer & mask))
				return true;
			mer >>= 2;
		}
		return is_poly(mer & mask); //last one
	}
};
#endif
