#ifndef _HASH_MAP_H
#define _HASH_MAP_H

#include <cstdint>
#ifdef _WIN32
#include <xmmintrin.h>
//#include <mmintrin.h>
#endif
#include <cstddef>

#include <algorithm>
#include <vector>
#include <utility>
#include <iterator>

#include "refresh/common/lib/defs.h"

namespace refresh
{
// ************************************************************************
// *** Global const iterator
template<typename HashMapLP>
class const_hash_map_lp_iterator
{
	friend HashMapLP;

public:
	using key_type = typename HashMapLP::key_type;
	using value_type = typename HashMapLP::value_type;
	using mapped_type = typename HashMapLP::mapped_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename HashMapLP::VectorType::const_iterator;
	using key_equal = typename HashMapLP::key_equal;

	const_hash_map_lp_iterator() = default;

	const_hash_map_lp_iterator(HashMapLP* _p_hm)
	{
		empty_key = _p_hm->empty_key;
		iter = _p_hm->data.end();
		iter_end = _p_hm->data.end();
		compare = _p_hm->compare;
	}

	const_hash_map_lp_iterator(HashMapLP* _p_hm, vector_iterator_type _iter)
	{
		iter = _iter;
		iter_end = _p_hm->data.end();
		empty_key = _p_hm->empty_key;
		compare = _p_hm->compare;
	}

	const value_type& operator*() const
	{
		return *iter;
	}

	const value_type* operator->() const
	{
		return &(*iter);
	}

	const_hash_map_lp_iterator<HashMapLP>& operator++()
	{
		increment();
		return *this;
	}

	const_hash_map_lp_iterator<HashMapLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	bool operator==(const const_hash_map_lp_iterator<HashMapLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	bool operator!=(const const_hash_map_lp_iterator<HashMapLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	vector_iterator_type iter;
	vector_iterator_type iter_end;
	key_type empty_key;
	key_equal compare;

	void increment()
	{
		for (++iter; iter != iter_end && compare(iter->first, empty_key); ++iter)
			;
	}
};

// ************************************************************************
// *** Global iterator
template <typename HashMapLP>
class hash_map_lp_iterator : public const_hash_map_lp_iterator<HashMapLP>
{
	friend HashMapLP;

public:
	using key_type = typename const_hash_map_lp_iterator<HashMapLP>::key_type;
	using value_type = typename const_hash_map_lp_iterator<HashMapLP>::value_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename HashMapLP::VectorType::iterator;

	hash_map_lp_iterator() = default;
	hash_map_lp_iterator(HashMapLP* _p_hm) : const_hash_map_lp_iterator<HashMapLP>(_p_hm)
	{
	}

	hash_map_lp_iterator(HashMapLP* _p_hm, vector_iterator_type _iter) : const_hash_map_lp_iterator<HashMapLP>(_p_hm, _iter)
	{
	}

	value_type& operator*()
	{
		return const_cast<value_type&>(*this->iter);
	}

	value_type* operator->()
	{
		return const_cast<value_type*>(&(*this->iter));
	}
};

// ************************************************************************
// *** Local const iterator
template<typename HashMapLP>
class const_hash_map_lp_local_iterator
{
	friend HashMapLP;

public:
	using key_type = typename HashMapLP::key_type;
	using value_type = typename HashMapLP::value_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename HashMapLP::VectorType::const_iterator;
	using vector_const_iterator_type = typename HashMapLP::VectorType::const_iterator;
	using key_equal = typename HashMapLP::key_equal;

	const_hash_map_lp_local_iterator() = default;

	const_hash_map_lp_local_iterator(HashMapLP* _p_hm)
	{
		//		local_key = _p_hm->empty_key;
		//		empty_key = _p_hm->empty_key;
		iter = _p_hm->data.end();
		//		iter_begin = _p_hm->data.begin();
		//		iter_end = _p_hm->data.end();
		//		compare = _p_hm->compare;
	}

	const_hash_map_lp_local_iterator(HashMapLP* _p_hm, vector_iterator_type _iter)
	{
		iter = _iter;
		iter_begin = _p_hm->data.begin();
		iter_end = _p_hm->data.end();
		empty_key = _p_hm->empty_key;
		if (iter != iter_end)
			local_key = iter->first;
		compare = _p_hm->compare;
	}

	const value_type& operator*() const
	{
		return *iter;
	}

	const value_type* operator->() const
	{
		return &(*iter);
	}

	const_hash_map_lp_local_iterator<HashMapLP>& operator++()
	{
		increment();
		return *this;
	}

	const_hash_map_lp_local_iterator<HashMapLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	bool operator==(const const_hash_map_lp_local_iterator<HashMapLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	bool operator!=(const const_hash_map_lp_local_iterator<HashMapLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	vector_iterator_type iter;
	vector_iterator_type iter_begin;
	vector_iterator_type iter_end;
	key_type local_key;
	key_type empty_key;
	key_equal compare;

	void increment()
	{
		for (++iter; iter != iter_end; ++iter)
		{
			if (compare(iter->first, empty_key))
			{
				iter = iter_end;
				return;
			}
			else if (compare(iter->first, local_key))
				return;
		}

		for (iter = iter_begin; iter != iter_end; ++iter)
		{
			if (compare(iter->first, empty_key))
			{
				iter = iter_end;
				return;
			}
			else if (compare(iter->first, local_key))
				return;
		}

		return;		// Should never be here
	}
};

// ************************************************************************
// *** Local const iterator
template<typename HashMapLP>
class const_hash_map_lp_local_value_iterator
{
	friend HashMapLP;

public:
	using key_type = typename HashMapLP::key_type;
	using value_type = typename HashMapLP::mapped_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename HashMapLP::VectorType::const_iterator;
	using vector_const_iterator_type = typename HashMapLP::VectorType::const_iterator;
	using key_equal = typename HashMapLP::key_equal;

	const_hash_map_lp_local_value_iterator() = default;

	const_hash_map_lp_local_value_iterator(HashMapLP* _p_hm)
	{
		//		local_key = _p_hm->empty_key;
		//		empty_key = _p_hm->empty_key;
		iter = _p_hm->data.end();
		//		iter_begin = _p_hm->data.begin();
		//		iter_end = _p_hm->data.end();
		//		compare = _p_hm->compare;
	}

	const_hash_map_lp_local_value_iterator(HashMapLP* _p_hm, vector_iterator_type _iter)
	{
		iter = _iter;
		iter_begin = _p_hm->data.begin();
		iter_end = _p_hm->data.end();
		empty_key = _p_hm->empty_key;
		if (iter != iter_end)
			local_key = iter->first;
		compare = _p_hm->compare;
	}

	const value_type& operator*() const
	{
		return iter->second;
	}

	const value_type* operator->() const
	{
		return &(iter->second);
	}

	const_hash_map_lp_local_value_iterator<HashMapLP>& operator++()
	{
		increment();
		return *this;
	}

	const_hash_map_lp_local_value_iterator<HashMapLP> operator++(int)
	{
		auto old_iter = *this;
		increment();
		return old_iter;
	}

	bool operator==(const const_hash_map_lp_local_value_iterator<HashMapLP>& rhs) const
	{
		return iter == rhs.iter;
	}

	bool operator!=(const const_hash_map_lp_local_value_iterator<HashMapLP>& rhs) const
	{
		return !(*this == rhs);
	}

protected:
	vector_iterator_type iter;
	vector_iterator_type iter_begin;
	vector_iterator_type iter_end;
	key_type local_key;
	key_type empty_key;
	key_equal compare;

	void increment()
	{
		for (++iter; iter != iter_end; ++iter)
		{
			if (compare(iter->first, empty_key))
			{
				iter = iter_end;
				return;
			}
			else if (compare(iter->first, local_key))
				return;
		}

		for (iter = iter_begin; iter != iter_end; ++iter)
		{
			if (compare(iter->first, empty_key))
			{
				iter = iter_end;
				return;
			}
			else if (compare(iter->first, local_key))
				return;
		}

		return;		// Should never be here
	}
};


// ************************************************************************
// *** Local iterator
template <typename HashMapLP>
class hash_map_lp_local_iterator : public const_hash_map_lp_local_iterator<HashMapLP>
{
	friend HashMapLP;

public:
	using key_type = typename const_hash_map_lp_iterator<HashMapLP>::key_type;
	using value_type = typename const_hash_map_lp_iterator<HashMapLP>::value_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename HashMapLP::VectorType::iterator;
	using vector_const_iterator_type = typename HashMapLP::VectorType::const_iterator;

	hash_map_lp_local_iterator() = default;
	hash_map_lp_local_iterator(HashMapLP* _p_hm) : const_hash_map_lp_local_iterator<HashMapLP>(_p_hm)
	{
	}

	hash_map_lp_local_iterator(HashMapLP* _p_hm, vector_iterator_type _iter) : const_hash_map_lp_local_iterator<HashMapLP>(_p_hm, _iter)
	{
	}

	hash_map_lp_local_iterator(HashMapLP* _p_hm, vector_const_iterator_type _iter) : const_hash_map_lp_local_iterator<HashMapLP>(_p_hm, _iter)
	{
	}

	value_type& operator*()
	{
		return const_cast<value_type&>(*this->iter);
	}

	value_type* operator->()
	{
		return const_cast<value_type*>(&(*this->iter));
	}
};

// ************************************************************************
// *** Local iterator
template <typename HashMapLP>
class hash_map_lp_local_value_iterator : public const_hash_map_lp_local_value_iterator<HashMapLP>
{
	friend HashMapLP;

public:
	using key_type = typename const_hash_map_lp_iterator<HashMapLP>::key_type;
	using value_type = typename const_hash_map_lp_iterator<HashMapLP>::mapped_type;
	using difference_type = ptrdiff_t;
	using iterator_category = std::forward_iterator_tag;
	using pointer = value_type*;
	using reference = value_type&;
	using vector_iterator_type = typename HashMapLP::VectorType::iterator;
	using vector_const_iterator_type = typename HashMapLP::VectorType::const_iterator;

	hash_map_lp_local_value_iterator() = default;
	hash_map_lp_local_value_iterator(HashMapLP* _p_hm) : const_hash_map_lp_local_value_iterator<HashMapLP>(_p_hm)
	{
	}

	hash_map_lp_local_value_iterator(HashMapLP* _p_hm, vector_iterator_type _iter) : const_hash_map_lp_local_value_iterator<HashMapLP>(_p_hm, _iter)
	{
	}

	hash_map_lp_local_value_iterator(HashMapLP* _p_hm, vector_const_iterator_type _iter) : const_hash_map_lp_local_value_iterator<HashMapLP>(_p_hm, _iter)
	{
	}

	value_type& operator*()
	{
		return const_cast<value_type&>(this->iter->second);
	}

	value_type* operator->()
	{
		return const_cast<value_type*>(&(this->iter->second));
	}
};

// ************************************************************************
// *** Hash map with linear probing (multikey)
template<typename Key_t, typename Value_t,
	typename Compare_t = std::equal_to<>,
	typename Hash_t = std::hash<Key_t>>
	class hash_map_lp {
	public:
		using key_type = Key_t;
		using mapped_type = Value_t;
		using value_type = std::pair<Key_t, Value_t>;
		using hasher = Hash_t;
		using key_equal = Compare_t;
		using reference = std::pair<Key_t, Value_t>&;
		using const_reference = const std::pair<Key_t, Value_t>&;
		using size_type = size_t;
		using difference_type = ptrdiff_t;
		using hash_map_lp_type = hash_map_lp<Key_t, Value_t, Compare_t, Hash_t>;
		using iterator = hash_map_lp_iterator<hash_map_lp_type>;
		using const_iterator = const_hash_map_lp_iterator<hash_map_lp_type>;
		using local_iterator = hash_map_lp_local_iterator<hash_map_lp_type>;
		using const_local_iterator = const_hash_map_lp_local_iterator<hash_map_lp_type>;
		using local_value_iterator = hash_map_lp_local_value_iterator<hash_map_lp_type>;
		using const_local_value_iterator = const_hash_map_lp_local_value_iterator<hash_map_lp_type>;

	private:
		using VectorType = std::vector<value_type>;

	private:
		Key_t empty_key;

		Compare_t compare;
		Hash_t hash;

		double max_fill_factor;
		size_t no_elements;
		size_t no_elements_unique;
		std::vector<value_type> data;
		size_t allocated;
		size_t size_when_restruct;
		size_t allocated_mask;

	public:
		friend class hash_map_lp_iterator<hash_map_lp_type>;
		friend class const_hash_map_lp_iterator<hash_map_lp_type>;
		friend class hash_map_lp_local_iterator<hash_map_lp_type>;
		friend class const_hash_map_lp_local_iterator<hash_map_lp_type>;
		friend class hash_map_lp_local_value_iterator<hash_map_lp_type>;
		friend class const_hash_map_lp_local_value_iterator<hash_map_lp_type>;

		virtual ~hash_map_lp() = default;

		explicit hash_map_lp(const Key_t _empty_key = Key_t(), size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Compare_t& _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			empty_key(_empty_key),
			compare(_compare),
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);
		}

		hash_map_lp(const hash_map_lp<Key_t, Value_t, Compare_t, Hash_t>& src) = default;

		hash_map_lp(hash_map_lp<Key_t, Value_t, Compare_t, Hash_t>&& src) = default;

		template <typename InputIterator>
		hash_map_lp(InputIterator first, InputIterator last,
			const Key_t _empty_key = Key_t(), size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Compare_t& _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			empty_key(_empty_key),
			compare(_compare),
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);

			for (auto p = first; p != last; ++p)
				insert_fast(*p);
		}

		explicit hash_map_lp(std::initializer_list<value_type> il,
			const Key_t _empty_key = Key_t(), size_t _init_reserved = 16, double _max_fill_factor = 0.7,
			const Compare_t& _compare = Compare_t(), const Hash_t _hash = Hash_t()) :
			empty_key(_empty_key),
			compare(_compare),
			hash(_hash)
		{
			construct(_init_reserved, _max_fill_factor);

			for (auto& x : il)
				insert_fast(x);
		}


		hash_map_lp_type& operator=(
			const hash_map_lp<Key_t, Value_t, Compare_t, Hash_t>& rhs)
		{
			if (this != &rhs)
			{
				data.clear();
				data.shrink_to_fit();

				compare = rhs.compare;
				hash = rhs.hash;
				data = rhs.data;
				max_fill_factor = rhs.max_fill_factor;
				no_elements = rhs.no_elements;
				no_elements_unique = rhs.no_elements_unique;
				allocated = rhs.allocated;
				size_when_restruct = rhs.size_when_restruct;
				allocated_mask = rhs.allocated_mask;
				empty_key = rhs.empty_key;
			}

			return *this;
		}

		hash_map_lp_type& operator=(
			hash_map_lp<Key_t, Value_t, Compare_t, Hash_t>&& rhs)
		{
			if (this != &rhs)
			{
				compare = rhs.compare;
				hash = rhs.hash;
				data = move(rhs.data);
				max_fill_factor = rhs.max_fill_factor;
				no_elements = rhs.no_elements;
				no_elements_unique = rhs.no_elements_unique;
				allocated = rhs.allocated;
				size_when_restruct = rhs.size_when_restruct;
				allocated_mask = rhs.allocated_mask;
				empty_key = rhs.empty_key;
			}

			return *this;
		}

		//	hash_map_type& operator=(std::initializer_list<value_type> il);

		void swap(hash_map_lp& x)
		{
			swap(compare, x.compare);
			swap(hash, x.hash);

			swap(max_fill_factor, x.max_fill_factor);
			swap(no_elements, x.no_elements);
			swap(no_elements_unique, x.no_elements_unique);
			swap(data, x.data);
			swap(allocated, x.allocated);
			swap(size_when_restruct, x.size_when_restruct);
			swap(allocated_mask, x.allocated_mask);

			swap(empty_key, x.empty_key);
		}

		// *************************************
		iterator begin()
		{
			if (!no_elements)
				return end();

			auto p = hash_map_lp_iterator<hash_map_lp_type>(this, data.begin());
			if (p->first == empty_key)
				++p;

			return p;
		}

		iterator end()
		{
			return hash_map_lp_iterator<hash_map_lp_type>(this);
		}

		const_iterator begin() const
		{
			return const_cast<hash_map_lp_type*>(this)->begin();
		}

		const_iterator end() const
		{
			return const_cast<hash_map_lp_type*>(this)->end();
		}

		const_iterator cbegin() const
		{
			return begin();
		}

		const_iterator cend() const
		{
			return end();
		}

		// *************************************
		local_iterator local_end()
		{
			return hash_map_lp_local_iterator<hash_map_lp_type>(this);
		}

		const_local_iterator local_end() const
		{
			return const_cast<hash_map_lp_type*>(this)->local_end();
		}

		/*		const_local_iterator clocal_end() const
				{
					return local_end();
				}*/

		local_iterator local_begin(iterator x)
		{
			return hash_map_lp_local_iterator<hash_map_lp_type>(this, x.iter);
		}

		local_iterator local_begin(const_iterator x)
		{
			return hash_map_lp_local_iterator<hash_map_lp_type>(this, x.iter);
		}

		// *************************************
		local_value_iterator local_value_end()
		{
			return hash_map_lp_local_value_iterator<hash_map_lp_type>(this);
		}

		const_local_value_iterator local_value_end() const
		{
			return const_cast<hash_map_lp_type*>(this)->local_value_end();
		}

		/*		const_local_iterator clocal_end() const
				{
					return local_end();
				}*/

		local_value_iterator local_value_begin(iterator x)
		{
			return hash_map_lp_local_value_iterator<hash_map_lp_type>(this, x.iter);
		}

		local_value_iterator local_value_begin(const_iterator x)
		{
			return hash_map_lp_local_value_iterator<hash_map_lp_type>(this, x.iter);
		}

		// *************************************
		bool empty()
		{
			return no_elements == 0;
		}

		size_type allocated_size()
		{
			return allocated;
		}

		size_type size() const
		{
			return no_elements;
		}

		size_type size_unique() const
		{
			return no_elements_unique;
		}

		size_type max_size()
		{
			return data.max_size();
		}

		size_t reserve(size_t _requested_reserve)
		{
			if (_requested_reserve < allocated)
				return allocated;

			restruct(_requested_reserve);
			return allocated;
		}

		//	T& operator[](const key_type& k);

		template <typename InputIterator>
		void insert(InputIterator first, InputIterator last)
		{
			for (auto p = first; p != last; ++p)
				insert_fast(*p);
		}

		void insert(std::initializer_list<value_type> il)
		{
			for (auto& x : il)
				insert_fast(x);
		}

		std::pair<iterator, bool> insert(const value_type& x)
		{
			if (no_elements >= size_when_restruct)
				restruct(allocated * 2);

			size_t h = hash(x.first) & allocated_mask;
			bool key_observed = false;

			if (!compare(data[h].first, empty_key))
			{
				key_observed |= compare(data[h].first, x.first);

				do
				{
					h = (h + 1) & allocated_mask;
					key_observed |= compare(data[h].first, x.first);
				} while (data[h].first != empty_key);
			}

			++no_elements;

			if (!key_observed)
				++no_elements_unique;

			data[h] = x;

			return std::make_pair(iterator(this, data.begin() + h), true);
		}

		bool insert_fast(const value_type& x)
		{
			if (no_elements >= size_when_restruct)
				restruct(allocated * 2);

			size_t h = hash(x.first) & allocated_mask;
			bool key_observed = false;

			if (!compare(data[h].first, empty_key))
			{
				key_observed |= compare(data[h].first, x.first);

				do
				{
					h = (h + 1) & allocated_mask;
					key_observed |= compare(data[h].first, x.first);
				} while (data[h].first != empty_key);
			}

			++no_elements;

			if (!key_observed)
				++no_elements_unique;

			data[h] = x;

			return true;
		}

		local_iterator find(const key_type& key)
		{
			// !!! For unknown reasons turning off inline here speeds things up significantly, but only in "find", it does not affect "check"
			size_t pos = _find_noinline(key, hash(key) & allocated_mask);

			if (compare(data[pos].first, key))
				return local_iterator(this, data.begin() + pos);
			else
				return local_iterator(this);
		}

		bool check(const key_type& key)
		{
			return _check_noinline(key, hash(key) & allocated_mask);

			/*			size_t pos = _find(key, hash(key) & allocated_mask);

						return compare(data[pos].first, key);*/
		}

		size_type count(const key_type& key) const
		{
			size_t pos = _find(key, hash(key) & allocated_mask);

			if (!compare(data[pos].key, key))
				return 0;

			size_t r = 0;

			for (; !compare(data[pos].key, empty_key); pos = (pos + 1) & allocated_mask)
				if (compare(data[pos].key, key))
					++r;

			return r;
		}

		key_equal key_eq() const
		{
			return compare;
		}

		hasher hash_function() const
		{
			return hash;
		}

		void prefetch(const Key_t& key)
		{
			size_t h = hash(key) & allocated_mask;

#ifdef _WIN32
			_mm_prefetch((const char*)(data.data() + h), _MM_HINT_T0);
#else
			__builtin_prefetch(data.data() + h);
#endif
		}

	private:
		void construct(size_t _init_reserved, double _max_fill_factor)
		{
			max_fill_factor = _max_fill_factor;

			if (max_fill_factor > 0.99)
				max_fill_factor = 0.99;
			else if (max_fill_factor < 0.01)
				max_fill_factor = 0.1;

			_reserve(_init_reserved);

			size_when_restruct = (size_t)(allocated * max_fill_factor);
			if (size_when_restruct == allocated)
				--size_when_restruct;
		}

		void _reserve(size_t requested_allocated)
		{
			allocated = requested_allocated;

			if (allocated < 8)
				allocated = 8;

			// Round up to the power of 2
			if ((allocated & (allocated - 1)))
			{
				while ((allocated & (allocated - 1)))
					allocated &= allocated - 1;
				allocated *= 2;
			}

			allocated_mask = allocated - 1ull;
			size_when_restruct = (size_t)(allocated * max_fill_factor);
			if (size_when_restruct == allocated)
				--size_when_restruct;

			value_type empty_cell;
			empty_cell.first = empty_key;

			data.resize(allocated, empty_cell);

			no_elements = 0;
			no_elements_unique = 0;
		}

		void restruct(size_t new_allocated)
		{
			std::vector<value_type> old_data;
			old_data = move(data);
			size_t old_allocated = allocated;

			_reserve(new_allocated);

			for (size_t i = 0; i < old_allocated; ++i)
				if (!compare(old_data[i].first, empty_key))
					insert_fast(old_data[i]);
		}

		size_t _find(const Key_t key, size_t pos)
		{
			if (compare(data[pos].first, empty_key) || compare(data[pos].first, key))
				return pos;

			for (++pos; pos < allocated; ++pos)
				if (compare(data[pos].first, empty_key) || compare(data[pos].first, key))
					return pos;

			for (pos = 0; pos < allocated; ++pos)
				if (compare(data[pos].first, empty_key) || compare(data[pos].first, key))
					return pos;

			return pos;
		}

		REFRESH_NO_INLINE size_t _find_noinline(const Key_t key, size_t pos)
		{
			if (compare(data[pos].first, empty_key) || compare(data[pos].first, key))
				return pos;

			for (++pos; pos < allocated; ++pos)
				if (compare(data[pos].first, empty_key) || compare(data[pos].first, key))
					return pos;

			for (pos = 0; pos < allocated; ++pos)
				if (compare(data[pos].first, empty_key) || compare(data[pos].first, key))
					return pos;

			return pos;
		}

		REFRESH_NO_INLINE bool _check_noinline(const Key_t key, size_t pos)
		{
			if (compare(data[pos].first, empty_key))
				return false;
			if (compare(data[pos].first, key))
				return true;

			for (++pos; pos < allocated; ++pos)
			{
				if (compare(data[pos].first, empty_key))
					return false;
				if (compare(data[pos].first, key))
					return true;
			}

			for (pos = 0; pos < allocated; ++pos)
			{
				if (compare(data[pos].first, empty_key))
					return false;
				if (compare(data[pos].first, key))
					return true;
			}

			return false;
		}

		bool _check(const Key_t key, size_t pos)
		{
			if (compare(data[pos].first, empty_key))
				return false;
			if (compare(data[pos].first, key))
				return true;

			for (++pos; pos < allocated; ++pos)
			{
				if (compare(data[pos].first, empty_key))
					return false;
				if (compare(data[pos].first, key))
					return true;
			}

			for (pos = 0; pos < allocated; ++pos)
			{
				if (compare(data[pos].first, empty_key))
					return false;
				if (compare(data[pos].first, key))
					return true;
			}

			return false;
		}
};

}
// EOF
#endif