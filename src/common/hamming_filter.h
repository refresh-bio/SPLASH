#ifndef _HAMMING_FILTER_H
#define _HAMMING_FILTER_H
#include <cinttypes>
#include <vector>
#include "target_count.h"
#include "../compactors/kmer.h"


class HammingFilter
{
    uint32_t min_hamming_threshold;

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
                if (static_cast<uint32_t>(KmerHelper::calculateHamming(target_i, target_j)) >= min_hamming_threshold)
                    return true;
            }
        }

        return false;
    }
};
#endif
