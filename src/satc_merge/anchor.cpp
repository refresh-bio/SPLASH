#include "anchor.h"
#include <cassert>
#include <algorithm>
#include <iostream>
#include <limits>
#include "../common/keep_n_largests.h"
#include "../common/binary_heap_merge.h"

//mkokot_TODO: consider remove from the code
Anchor merge_keep_target_order(const std::vector<Anchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt) {
	assert(to_merge.size());
	n_unique_targets = 0;
	tot_cnt = 0;
	Anchor res;
	res.anchor = to_merge[0].anchor;
	std::vector<uint64_t> read_pos(to_merge.size());

	size_t records_left = 0;
	for (const auto& x : to_merge)
		records_left += x.data.size();

	while (records_left) {	// !!! TODO: dla sporej liczby sampli pewnie lepsze by bylo scalanie z uzyciem kopca
		uint64_t min_target = std::numeric_limits<uint64_t>::max();

		//find min target
		for (size_t i = 0; i < to_merge.size(); ++i) {
			if (read_pos[i] < to_merge[i].data.size()) {
				if (to_merge[i].data[read_pos[i]].target < min_target) {
					min_target = to_merge[i].data[read_pos[i]].target;
				}
			}
		}

		++n_unique_targets;

		//store all record having min target in the result
		for (size_t i = 0; i < to_merge.size(); ++i) {
			if (read_pos[i] < to_merge[i].data.size()) {
				if (to_merge[i].data[read_pos[i]].target == min_target) {
					res.data.push_back(to_merge[i].data[read_pos[i]]);
					tot_cnt += res.data.back().count;

					++read_pos[i];
					--records_left;
				}
			}
		}
	}
	res.data.shrink_to_fit();
	return res;
}
//mkokot_TODO: consider remove from the code
Anchor merge_keep_target_order_binary_heap(const std::vector<Anchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt) {
	assert(to_merge.size());
	n_unique_targets = 0;
	tot_cnt = 0;
	Anchor res;
	res.anchor = to_merge[0].anchor;
	std::vector<uint64_t> read_pos(to_merge.size());

	struct heap_desc_t
	{
		uint64_t elem;
		uint64_t id;
		heap_desc_t(uint64_t elem, uint64_t id) :
			elem(elem), id(id) {

		}
		bool operator<(const heap_desc_t& rhs) const {
			return elem > rhs.elem;
		}
	};
	
	std::vector<heap_desc_t> heap;
	heap.reserve(to_merge.size());

	for (size_t id = 0; id < to_merge.size(); ++id) {
		if (to_merge[id].data.size())
			heap.emplace_back(to_merge[id].data[0].target, id);
	}

	std::make_heap(heap.begin(), heap.end());


	auto heap_down = [&heap]() {
		uint64_t parent = 0;
		uint64_t left = 1;
		uint64_t right = 2;


		auto elem = heap[0];

		//has 2 childs
		while (right < heap.size()) {
			left = heap[left].elem < heap[right].elem ? left : right; //left is min now

			if (elem.elem <= heap[left].elem) {
				heap[parent] = elem;
				return;
			}

			heap[parent] = heap[left];
			parent = left;
			left = parent * 2 + 1;
			right = left + 1;
		}

		if (left < heap.size()) {
			if (elem.elem < heap[left].elem) {
				heap[parent] = elem;
				return;
			}
			heap[parent] = heap[left];
			parent = left;
		}

		heap[parent] = elem;
	};

	//first element
	uint64_t target = heap[0].elem;

	uint64_t id = heap[0].id;
	++n_unique_targets;
	
	tot_cnt += to_merge[id].data[read_pos[id]].count;

	res.data.push_back(to_merge[id].data[read_pos[id]++]);

	if (read_pos[id] < to_merge[id].data.size())
	{
		heap[0].elem = to_merge[id].data[read_pos[id]].target;
		heap_down();
	}
	else
	{
		heap[0] = heap.back();
		heap.pop_back();
		if (heap.size() > 1)
			heap_down();
	}
	

	while (heap.size()) {
		auto new_target = heap[0].elem;		
		id = heap[0].id;

		tot_cnt += to_merge[id].data[read_pos[id]].count;

		res.data.push_back(to_merge[id].data[read_pos[id]++]);

		if (read_pos[id] < to_merge[id].data.size())
		{
			heap[0].elem = to_merge[id].data[read_pos[id]].target;
			heap_down();
		}
		else
		{
			heap[0] = heap.back();
			heap.pop_back();
			if (heap.size() > 1)
				heap_down();
		}
		

		if (new_target != target) {
			++n_unique_targets;
			target = new_target;
		}
	}
	return res;

}


Anchor merge_keep_target_order_binary_heap(const std::vector<Non10SingleSampleAnchor>& to_merge, uint64_t& n_unique_targets, uint64_t& tot_cnt) {
	assert(to_merge.size());
	n_unique_targets = 0;
	tot_cnt = 0;
	Anchor res;
	res.anchor = to_merge[0].anchor;

	auto get_elem = [&](size_t id, size_t pos) {
		return to_merge[id].data[pos].target;
	};
	auto get_array_size = [&](size_t id) {
		return to_merge[id].data.size();
	};

	BinaryHeapMerge<uint64_t> heap(to_merge.size(), get_elem, get_array_size);

	uint64_t cur_target;

	//first elem
	heap.ProcessElem(get_elem, get_array_size, [&](uint64_t target, size_t id, size_t pos) {
		cur_target = target;
		++n_unique_targets;
		tot_cnt += to_merge[id].data[pos].count;
		res.data.emplace_back(
			0,
			to_merge[id].data[pos].target,
			to_merge[id].sample_id,
			to_merge[id].data[pos].count
		);
	});

	while (!heap.Empty()) {

		heap.ProcessElem(get_elem, get_array_size, [&](uint64_t target, size_t id, size_t pos) {
			tot_cnt += to_merge[id].data[pos].count;

			res.data.emplace_back(
				0,
				to_merge[id].data[pos].target,
				to_merge[id].sample_id,
				to_merge[id].data[pos].count
			);

			if (target != cur_target) {
				++n_unique_targets;
				cur_target = target;
			}
		});
	}
	return res;
}

Anchor merge_keep_target_order_binary_heap(const std::vector<Non10SingleSampleAnchor>& to_merge, uint64_t keep_n_most_freq_targets, uint64_t& n_unique_targets, uint64_t& tot_cnt, uint64_t& n_unique_targets_kept, uint64_t& tot_cnt_kept) {
	assert(to_merge.size());
	if (keep_n_most_freq_targets == 0) {
		auto res = merge_keep_target_order_binary_heap(to_merge, n_unique_targets, tot_cnt);
		n_unique_targets_kept = n_unique_targets;
		tot_cnt_kept = tot_cnt;
		return res;
	}

	assert(keep_n_most_freq_targets);
	n_unique_targets = 0;
	tot_cnt = 0;
	Anchor res;
	res.anchor = to_merge[0].anchor;

	auto get_elem = [&](size_t id, size_t pos) {
		return to_merge[id].data[pos].target;
	};
	auto get_array_size = [&](size_t id) {
		return to_merge[id].data.size();
	};

	BinaryHeapMerge<uint64_t> heap(to_merge.size(), get_elem, get_array_size);

	struct elem_t {
		uint64_t target{};
		uint64_t cnt{};
		struct sample_id_cnt {
			uint64_t sample_id;
			uint64_t cnt;
			sample_id_cnt(uint64_t sample_id, uint64_t cnt) : sample_id(sample_id), cnt(cnt) {

			}
		};
		std::vector<sample_id_cnt> sample_ids_cnts;
		bool operator>(const elem_t& rhs) const {
			return cnt > rhs.cnt;
		}
	};
	elem_t e;

	KeepNLargests<elem_t> filtering_heap(keep_n_most_freq_targets);

	//first elem
	heap.ProcessElem(get_elem, get_array_size, [&](uint64_t target, size_t id, size_t pos) {
		++n_unique_targets;
		tot_cnt += to_merge[id].data[pos].count;

		e.target = target;
		e.cnt = to_merge[id].data[pos].count;
		e.sample_ids_cnts.emplace_back(to_merge[id].sample_id, to_merge[id].data[pos].count);
	});

	while (!heap.Empty()) {

		heap.ProcessElem(get_elem, get_array_size, [&](uint64_t target, size_t id, size_t pos) {
			tot_cnt += to_merge[id].data[pos].count;

			if (target != e.target) {
				filtering_heap.Add(std::move(e));
				++n_unique_targets;

				e.target = target;
				e.cnt = to_merge[id].data[pos].count;
				e.sample_ids_cnts.clear();
				e.sample_ids_cnts.emplace_back(to_merge[id].sample_id, to_merge[id].data[pos].count);
			} else {
				e.cnt += to_merge[id].data[pos].count;
				e.sample_ids_cnts.emplace_back(to_merge[id].sample_id, to_merge[id].data[pos].count);
			}
		});
	}

	filtering_heap.Add(std::move(e));

	std::vector<elem_t> filtered;
	filtering_heap.StealSorted(filtered, [](const elem_t& lhs, const elem_t& rhs) {
		return lhs.target < rhs.target;
	});

	//repack
	n_unique_targets_kept = filtered.size();
	tot_cnt_kept = 0;

	for (const elem_t& elem : filtered) {
		tot_cnt_kept += elem.cnt;
		for (auto sample_id_cnt : elem.sample_ids_cnts) {
			res.data.emplace_back(
				0,
				elem.target,
				sample_id_cnt.sample_id,
				sample_id_cnt.cnt
			);
		}
	}

	return res;
}

Anchor keep_n_most_freq_targets(const Anchor& input, uint64_t keep_n_most_freq_targets) {
	//this should not be called if keep_n_most_freq_targets is zero which means keep all
	assert(keep_n_most_freq_targets);
	//this also should be checked outside to prevent unnecessary copy
	assert(input.data.size() > keep_n_most_freq_targets);

	Anchor res;
	res.anchor = input.anchor;

	struct elem_t {
		uint64_t target{};
		uint64_t cnt{};
		struct sample_id_barcode_cnt {
			uint64_t sample_id;
			uint64_t barcode;
			uint64_t cnt;
			sample_id_barcode_cnt(uint64_t sample_id, uint64_t barcode, uint64_t cnt) :
				sample_id(sample_id),
				barcode(barcode),
				cnt(cnt) {

			}
		};
		std::vector<sample_id_barcode_cnt> sample_ids_barcodes_cnts;
		bool operator>(const elem_t& rhs) const {
			return cnt > rhs.cnt;
		}
	};
	elem_t e;

	KeepNLargests<elem_t> heap(keep_n_most_freq_targets);

	e.target = input.data[0].target;
	e.cnt = input.data[0].count;
	e.sample_ids_barcodes_cnts.emplace_back(input.data[0].sample_id, input.data[0].barcode, input.data[0].count);

	for (uint64_t i = 1 ; i < input.data.size(); ++i) {
		if (input.data[i].target != e.target) {
			//just checking if targets are sorted
			assert(input.data[i].target > e.target);

			heap.Add(std::move(e));

			e.target = input.data[i].target;
			e.cnt = input.data[i].count;
			e.sample_ids_barcodes_cnts.clear();
			e.sample_ids_barcodes_cnts.emplace_back(input.data[i].sample_id, input.data[i].barcode, input.data[i].count);
		} else {
			e.cnt += input.data[i].count;
			e.sample_ids_barcodes_cnts.emplace_back(input.data[i].sample_id, input.data[i].barcode, input.data[i].count);
		}
	}

	heap.Add(std::move(e));

	std::vector<elem_t> to_keep;
	heap.StealSorted(to_keep, [](const elem_t& lhs, const elem_t& rhs) {
		return lhs.target < rhs.target;
	});

	//repack
	for (elem_t& elem : to_keep) {
		//sorting because we store in Header::ordering_t::ATSBC
		std::sort(elem.sample_ids_barcodes_cnts.begin(), elem.sample_ids_barcodes_cnts.end(), 
			[](const elem_t::sample_id_barcode_cnt& lhs, const elem_t::sample_id_barcode_cnt& rhs)
			{
				if (lhs.sample_id != rhs.sample_id)
					return lhs.sample_id < rhs.sample_id;
				if (lhs.barcode != rhs.barcode)
					return lhs.barcode < rhs.barcode;
				return lhs.cnt < rhs.cnt;
			});

		for (auto sample_id_barcode_cnt : elem.sample_ids_barcodes_cnts) {
			res.data.emplace_back(
				sample_id_barcode_cnt.barcode,
				elem.target,
				sample_id_barcode_cnt.sample_id,
				sample_id_barcode_cnt.cnt
			);
		}
	}
	return res;
}
