#ifndef _TARGET_COUNT_H
#define _TARGET_COUNT_H
#include <cstdint>

struct TargetCount {
	uint64_t target;
	uint32_t count;
	TargetCount(uint64_t target, uint32_t count) :
		target(target),
		count(count)
	{

	}
};

#endif //_TARGET_COUNT_H
