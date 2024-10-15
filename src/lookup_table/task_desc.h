#ifndef _PACK_DESC
#define _PACK_DESC
#include <cstdint>

struct task_desc {

	uint64_t kmer{}; //always
	uint32_t file_id{}; // cat 2
	uint32_t header_id{}; //cat 1 or 3
	uint32_t cnt{}; // cat 2

	//result of querying kmc databases
	int cat; //1 or 2 or 3, -1 means unknown yet

	int64_t idx = -1; //result

	task_desc(uint64_t kmer, uint32_t file_id, uint32_t header_id, uint32_t cnt = 0, int cat = -1) :
		kmer(kmer), file_id(file_id), header_id(header_id), cnt(cnt), cat(cat) {
	}
};

#endif // _PACK_DESC
