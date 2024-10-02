#ifndef _BINARY_HEAP_MERGE_H
#define _BINARY_HEAP_MERGE_H

#include <vector>

template<typename T>
class BinaryHeapMerge {
	template<typename U>
	struct heap_desc_t
	{
		U elem;
		uint64_t id;
		heap_desc_t(U elem, uint64_t id) :
			elem(elem), id(id) {

		}
		bool operator<(const heap_desc_t& rhs) const {
			return elem > rhs.elem;
		}
	};

	std::vector<size_t> read_pos;
	std::vector<heap_desc_t<T>> heap;

	void heap_down() {
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

public:
	template<typename GET_ELEM, typename GET_ARRAY_SIZE>
	BinaryHeapMerge(size_t n_arrays_to_merge, GET_ELEM&& get_elem, GET_ARRAY_SIZE&& get_array_size) :
		read_pos(n_arrays_to_merge) {
		heap.reserve(n_arrays_to_merge);

		for (size_t id = 0; id < n_arrays_to_merge; ++id) {
			if (get_array_size(id))
				heap.emplace_back(get_elem(id, 0), id);
		}

		std::make_heap(heap.begin(), heap.end());
	}

	template<typename GET_ELEM, typename GET_ARRAY_SIZE, typename Callback>
	void ProcessElem(GET_ELEM&& get_elem, GET_ARRAY_SIZE&& get_array_size, Callback&& callback) {
		auto id = heap[0].id;

		callback(heap[0].elem, id, read_pos[id]);

		if (++read_pos[id] < get_array_size(id))
		{
			heap[0].elem = get_elem(id, read_pos[id]);
			heap_down();
		}
		else
		{
			heap[0] = heap.back();
			heap.pop_back();
			if (heap.size() > 1)
				heap_down();
		}
	}
	bool Empty() const {
		return heap.empty();
	}
};



//when we don't know the total size
//in fact the above could be implemented using this as implementation and just wrapping it with vector of read_pos
template<typename T>
class BinaryHeapMergeStreams {
	template<typename U>
	struct heap_desc_t
	{
		U elem;
		uint64_t id;
		heap_desc_t(U elem, uint64_t id) :
			elem(elem), id(id) {

		}
		bool operator<(const heap_desc_t& rhs) const {
			return elem > rhs.elem;
		}
	};

	std::vector<heap_desc_t<T>> heap;

	void heap_down() {
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

public:
	template<typename DO_WITH_ELEM_IF_EXISTS>
	BinaryHeapMergeStreams(size_t n_streams_to_merge, const DO_WITH_ELEM_IF_EXISTS& do_with_elem_if_exists) {
		heap.reserve(n_streams_to_merge);

		for (size_t id = 0; id < n_streams_to_merge; ++id) {
			do_with_elem_if_exists(id, [&](const T& elem) { heap.emplace_back(elem, id); });
		}

		std::make_heap(heap.begin(), heap.end());
	}

	template<typename DO_WITH_ELEM_IF_EXISTS, typename Callback>
	void ProcessElem(const DO_WITH_ELEM_IF_EXISTS& do_with_elem_if_exists, Callback&& callback) {
		auto id = heap[0].id;

		callback(heap[0].elem, id);

		if (do_with_elem_if_exists(id, [&](const T& elem)
			{
				heap[0].elem = elem;
				heap_down();
			}))
		{} //do nothing in if's body, because it is done in the callback of do_with_elem_if_exists
		else
		{
			heap[0] = heap.back();
			heap.pop_back();
			if (heap.size() > 1)
				heap_down();
		}
	}
	bool Empty() const {
		return heap.empty();
	}
};

#endif // !_BINARY_HEAP_MERGE_H