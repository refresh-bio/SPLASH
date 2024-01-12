#ifndef _KEEP_N_LARGESTS
#define _KEEP_N_LARGESTS
#include <vector>
#include <algorithm>

template<typename T, typename PRED = std::greater<T>>
class KeepNLargests {
	std::vector<T> heap;
	size_t n;
	PRED pred;	
public:
	KeepNLargests(size_t n, PRED pred = PRED{}) : n(n), pred(pred) {

	}
	void Add(T&& elem) {
		if (heap.size() == n) {
			if (pred(heap[0], elem))
				return;
			heap.push_back(std::move(elem));
			std::push_heap(heap.begin(), heap.end(), pred);
			std::pop_heap(heap.begin(), heap.end(), pred);
			heap.pop_back();
		}
		else {
			heap.push_back(std::move(elem));
			if (heap.size() == n) {
				std::make_heap(heap.begin(), heap.end(), pred);
			}
		}
	}
	const std::vector<T>& Get() const {
		return heap;
	}

	void Steal(std::vector<T>& res) const {
		res = std::move(heap);
	}

	template<typename Pred = std::less<T>>
	std::vector<T> GetSorted(Pred pred = std::less<T>{}) {
		std::vector<T> res = heap;
		std::sort(res.begin(), res.end(), pred);
		return res;
	}

	template<typename Pred = std::less<T>>
	void StealSorted(std::vector<T>& res, Pred pred = std::less<T>{}) {
		std::sort(heap.begin(), heap.end(), pred);
		res = std::move(heap);
	}
};

#endif // !_KEEP_N_LARGESTS
