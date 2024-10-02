#pragma once
#include <refresh/memory_chunk/lib/memory_chunk.h>

//mkokot_TODO: move it to devlibs? consider name and if we want an interface to extend capacity
// meaning we will enlarge keeping current data (up to size or up to capacity?)
// or maybe we should behave just like vector, but also add unsafa methods upush_back (unsafe push back that does not check size < capacity and reallocate), uresize (unsafe resize), etc.
//const capacity vector with no checks
template<typename T>
class unsafe_vector {
	std::unique_ptr<T[]> _raw_data;
	refresh::memory_chunk<T> chunk;
public:
	using size_type = std::size_t;
	using value_type = T;
	using reference = value_type&;
	using const_reference = const value_type&;
	using iterator = value_type*;
	using const_iterator = const value_type*;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;
	using pointer = value_type*;
	using const_pointer = const value_type*;

	unsafe_vector() = default;

	unsafe_vector(size_t capacity) :
		_raw_data(new T[capacity]),   //_raw_data(std::make_unique_for_overwrite<T[]>(capacity)), //in c++20
		chunk(_raw_data.get(), capacity) {
	}

	void clear() {
		chunk.clear();
	}
	void resize(size_type new_size) {
		chunk.resize(new_size);
	}
	void reserve(size_type new_capacity) {
		chunk.reserve(new_capacity);
	}
	size_type size() const {
		return chunk.size();
	}
	size_type capacity() const {
		return chunk.capacity();
	}
	size_type max_size() const {
		return chunk.max_size();
	}
	reference operator[](size_type pos) {
		return chunk[pos];
	}
	const_reference operator[](size_type pos) const {
		return chunk[pos];
	}
	reference at(size_type pos)
	{
		return chunk.at(pos);
	}
	const_reference at(size_type pos) const
	{
		return chunk.at(pos);
	}
	void push_back(const T& val)
	{
		chunk.push_back(val);
	}
	void push_back(T&& val)
	{
		chunk.emplace_back(std::move(val));
	}
	template< class... Args >
	reference emplace_back(Args&&... args)
	{
		chunk.emplace_back(std::forward<Args>(args)...);
	}
	void pop_back()
	{
		chunk.pop_back();
	}
	bool empty() const
	{
		return chunk.empty();
	}
	reference front()
	{
		return chunk.front();
	}
	reference back()
	{
		return chunk.back();
	}
	const_reference front() const
	{
		return chunk.front();
	}
	const_reference back() const
	{
		return chunk.back();
	}
	iterator begin()
	{
		return chunk.begin();
	}
	iterator end()
	{
		return chunk.end();
	}
	const_iterator begin() const
	{
		return chunk.begin();
	}
	const_iterator end() const
	{
		return chunk.end();
	}
	const_iterator cbegin()
	{
		return chunk.cbegin();
	}
	const_iterator cend()
	{
		return chunk.cend();
	}
	const_iterator cbegin() const
	{
		return chunk.cbegin();
	}
	const_iterator cend() const
	{
		return chunk.cend();
	}
	reverse_iterator rbegin()
	{
		return chunk.rbegin();
	}
	reverse_iterator rend()
	{
		return chunk.rend();
	}
	const_reverse_iterator rbegin() const
	{
		return chunk.rbegin();
	}
	const_reverse_iterator rend() const
	{
		return chunk.rend();
	}
	const_reverse_iterator crbegin()
	{
		return chunk.crbegin();
	}
	const_reverse_iterator crend()
	{
		return chunk.crend();
	}
	const_reverse_iterator crbegin() const
	{
		return chunk.crbegin();
	}
	const_reverse_iterator crend() const
	{
		return chunk.crend();
	}
	pointer data()
	{
		return chunk.data();
	}
	const_pointer data() const
	{
		return chunk.data();
	}
};
