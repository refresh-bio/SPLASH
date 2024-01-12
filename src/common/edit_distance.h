#pragma once
#include <cstdint>

class CEditDistanceOneWord
{
	uint32_t len;
	uint64_t masks[4];

public:
	CEditDistanceOneWord() = default;

	void Prepare(const uint64_t seq1, uint32_t _len);
	uint32_t Calculate(const uint64_t seq2);
	uint32_t Calculate(const uint64_t seq1, const uint64_t seq2, uint32_t _len);
};