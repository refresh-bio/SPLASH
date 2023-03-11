#include "edit_distance.h"


// *********************************************************************************************
void CEditDistanceOneWord::Prepare(const uint64_t seq1, uint32_t _len)
{
	uint64_t x = seq1;
	len = _len;

	masks[0] = masks[1] = masks[2] = masks[3] = 0;

	for (uint32_t i = 0; i < len; ++i)
	{
		masks[x & 3ull] += 1ull << i;
		x >>= 2;
	}
}

// *********************************************************************************************
// Based on Hyyro's bit-par algorithm "A Bit-Vector Algorithm for Computing Levenshtein and Damerau Edit Distances"
// https://www.semanticscholar.org/paper/A-Bit-Vector-Algorithm-for-Computing-Levenshtein-Hyyr%C3%B6/813e26d8920d17c2afac6bf5a15c537b067a128a/figure/1
uint32_t CEditDistanceOneWord::Calculate(const uint64_t seq2)
{
	uint64_t VN = 0;
	uint64_t VP = ~0ull >> (64 - len);
	uint64_t test_bit = 1ull << (len - 1);

	uint64_t Rp, HPp, HNp, VPp, VNp;
	uint32_t ed = len;

	uint64_t x = seq2;

	for (uint32_t i = 0; i < len; ++i)
	{
		uint64_t symbol = x & 3ull;
		x >>= 2;

		Rp = (((masks[symbol] & VP) + VP) ^ VP) | masks[symbol] | VN;
		HPp = VN | ~(Rp | VP);
		HNp = VP & Rp;
		//		VPp = (HNp << 1) | ~(Rp | (HPp << 1));
		//		VNp = (HPp << 1) & Rp;
		VPp = (HNp << 1) | ~(Rp | ((HPp << 1) | 1));
		VNp = ((HPp << 1) | 1) & Rp;

		if (HPp & test_bit)
			ed++;
		else if (HNp & test_bit)
			ed--;

		VP = VPp;
		VN = VNp;
	}

	return ed;
}

// *********************************************************************************************
uint32_t CEditDistanceOneWord::Calculate(const uint64_t seq1, const uint64_t seq2, uint32_t _len)
{
	Prepare(seq1, _len);

	return Calculate(seq2);
}