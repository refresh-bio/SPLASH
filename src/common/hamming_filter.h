#ifndef _HAMMING_FILTER_H
#define _HAMMING_FILTER_H
#include <cinttypes>
#include <vector>
#include "target_count.h"

class HammingFilter
{
    uint32_t min_hamming_threshold;

 static int hamming_dist(uint64_t a, uint64_t b) {
		// generate mask with even bits
		static uint64_t MASK_EVEN_BITS = []() {
			uint64_t v = 0;
			for (int i = 0; i < 32; ++i) {
				v <<= 2;
				v |= 1ULL;
			}
			return v;
		}();

		uint64_t diffs_even = (a ^ b) & MASK_EVEN_BITS;
		uint64_t diffs_odd = ((a >> 1) ^ (b >> 1)) & MASK_EVEN_BITS;

		uint64_t diffs = diffs_even | diffs_odd;
		#ifdef _WIN32
		int cnt = __popcnt64(diffs);
		#else
		int cnt = __builtin_popcountll(diffs);
		#endif
		return cnt;
	}

public:
    HammingFilter(uint32_t min_hamming_threshold)
        : min_hamming_threshold(min_hamming_threshold)
    {

    }
    HammingFilter() :
        HammingFilter(0)
    {

    }

    bool ContainsDistantPair(const std::vector<TargetCount>& targets) const
    {
        if (min_hamming_threshold == 0)
            return true;

        if (targets.size() < 2)
            return false;

        for (size_t i = 0; i < targets.size(); ++i)
        {
            for (size_t j = i+1; j < targets.size(); ++j)
            {
                uint64_t target_i = targets[i].target;
                uint64_t target_j = targets[j].target;
                if (static_cast<uint32_t>(hamming_dist(target_i, target_j)) >= min_hamming_threshold)
                    return true;
            }
        }

        return false;
    }
};
#endif
