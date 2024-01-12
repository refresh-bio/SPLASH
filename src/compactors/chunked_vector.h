#pragma once
#include <vector>
#include <deque>



template <class T>
class chunked_vector {

	std::vector<std::vector<T>> chunks;
	size_t chunk_size{ 0 };
	size_t chunk_id_bits{ 0 };
	size_t chunk_offset_mask{ 0 };

	const std::vector<T>& last_chunk() const { return chunks.back(); }

	std::vector<T>& last_chunk() { return chunks.back(); }

	void extend() {
		chunks.push_back(std::vector<T>());
		chunks.back().reserve(chunk_size);
	}

public:
	chunked_vector(size_t chunk_size) {

		// round chunk size up to the closest power of two and fill masks
		while (chunk_size > 0) {
			++chunk_id_bits;
			chunk_offset_mask <<= 1;
			chunk_offset_mask |= 1;
			chunk_size >>= 1;
		}

		this->chunk_size = chunk_offset_mask + 1;
		extend();
	}

	~chunked_vector() {
	}

	size_t num_chunks() const { return chunks.size(); }

	size_t capacity() const {
		return chunks.size() * chunk_size;
	}

	T* resize_for_additional(size_t num) {
		size_t new_size = last_chunk().size() + num;
		
		if (new_size >= chunk_size) {
			chunks.push_back(std::vector<T>());
			new_size = num;
			last_chunk().reserve(std::max(chunk_size, new_size));
		}

		T* place = last_chunk().data() + last_chunk().size();
		last_chunk().resize(new_size);

		return place;
	}

	template<class... Args>
	void emplace_back(Args... args) {

		last_chunk().emplace_back(args...);

		if (last_chunk().size() == chunk_size) {
			extend();
		}
	}

	const T& operator [](size_t i) const {
		return chunks[i >> chunk_id_bits][i & chunk_offset_mask];
	}

	T& operator [](size_t i) {
		return chunks[i >> chunk_id_bits][i & chunk_offset_mask];
	}

	T& back() { return  chunks.back().back(); }

	const T& back() const { return chunks.back().back(); }
};


/*
template <class T>
class flexible_deque {
	
	std::vector<T*> chunks;
	size_t chunk_size{ 0 };
	size_t chunk_id_bits{ 0 };
	size_t chunk_offset_mask{ 0 };

	size_t topOffset{ 0 };


public:
	flexible_deque(size_t chunk_size) {
	
		// round chunk size up to the closest power of two and fill masks
		while (chunk_size > 0) {
			++chunk_id_bits;
			chunk_offset_mask <<= 1;
			chunk_offset_mask |= 1;
			chunk_size >>= 1;
		}

		this->chunk_size = chunk_offset_mask + 1;
		reserve(this->chunk_size);
	}

	~flexible_deque() {
		for (auto c : chunks) {
			delete[] c;
		}
	}

	const T& operator [](size_t i) const {
		return chunks[i >> chunk_id_bits][i & chunk_offset_mask];
	}

	T& operator [](size_t i) {
		return chunks[i >> chunk_id_bits][i & chunk_offset_mask];
	}

	size_t capacity() const {
		return chunks.size() * chunk_size;
	}

	void reserve(size_t size) {
		while (capacity() < size) {
			chunks.push_back(new T[chunk_size]);
		}
	}

	T* push(const T& v) {
		
		T* place = &chunks.back()[topOffset];
		*place = v;

		if (++topOffset == chunk_size) {
			chunks.push_back(new T[chunk_size]);
			topOffset = 0;
		}

		return place;
	}

	T* push(const T* vals, int num) {
		// check if we go beyond (it may leave one empty slot at the very end)
		if (topOffset + num >= chunk_size) {
			chunks.push_back(new T[chunk_size]);
			topOffset = 0;
		}

		T* place = &chunks.back()[topOffset];
		
		for (int i = 0; i < num; ++i) {
			place[i] = vals[i];
		}

		topOffset += num;
		return place;
	}
};
*/