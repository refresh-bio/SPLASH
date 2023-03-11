#ifndef _ANCHOR_H
#define _ANCHOR_H

#include <vector>
#include <cinttypes>

struct AnchorData {
	uint64_t barcode;
	uint64_t target;
	uint64_t sample_id;
	uint64_t count;
	AnchorData(uint64_t barcode, uint64_t target, uint64_t sample_id, uint64_t count) :
		barcode(barcode),
		target(target),
		sample_id(sample_id),
		count(count)
	{

	}
};

struct Anchor {
	uint64_t anchor;

	std::vector<AnchorData> data;
};

struct Non10SingleSampleAnchorData {
	uint64_t target;
	uint64_t count;
	Non10SingleSampleAnchorData(uint64_t target, uint64_t count) :
		target(target),
		count(count)
	{

	}
};

struct Non10SingleSampleAnchor {
	uint64_t anchor;
	uint64_t sample_id;

	std::vector<Non10SingleSampleAnchorData> data;
};

Anchor merge_keep_target_order(const std::vector<Anchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt);

Anchor merge_keep_target_order_binary_heap(const std::vector<Anchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt);

Anchor merge_keep_target_order_binary_heap(const std::vector<Non10SingleSampleAnchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt);

Anchor merge_keep_target_order_binary_heap(const std::vector<Non10SingleSampleAnchor>& to_merge, uint64_t keep_n_most_freq_targets, uint64_t& n_unique_targets, uint64_t& tot_cnt, uint64_t& n_unique_targets_kept, uint64_t& tot_cnt_kept);

#endif // _ANCHOR_H

