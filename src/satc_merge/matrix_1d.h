#pragma once
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

namespace refresh {

template<typename T> class matrix_1d
{
	T* ptr;
	size_t n_elems;

	void allocate()
	{
		ptr = new T[n_elems];
	}

	void release()
	{
		if (ptr)
		{
			delete[] ptr;
			ptr = nullptr;

			n_elems = 0;
		}
	}

public:
	matrix_1d() : ptr(nullptr), n_elems(0) {}

	matrix_1d(const size_t n_elems) : ptr(nullptr), n_elems(n_elems)
	{
		allocate();
	}

	matrix_1d(const size_t n_elems, const T val) : ptr(nullptr), n_elems(n_elems)
	{
		allocate();
		std::fill_n(ptr, n_elems, val);
	}

	matrix_1d(const matrix_1d<T>& x)
	{
		n_elems = x.n_elems;

		allocate();

		std::memcpy(ptr, x.ptr, n_elems * sizeof(T));
	}

	matrix_1d(matrix_1d<T>&& x)
	{
		n_elems = x.n_elems;
		ptr = x.ptr;

		x.n_elems = 0;
		x.ptr = nullptr;
	}

	matrix_1d(const std::initializer_list<T>& il)
	{
		n_elems = il.size();
		ptr = nullptr;
		allocate();

		auto p = il.begin();

		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] = *p++;
	}

	~matrix_1d()
	{
		release();
	}

	template<typename Iter>
	matrix_1d(const Iter begin, const Iter end)
	{
		n_elems = std::distance(begin, end);
		ptr = nullptr;
		allocate();

		auto q = ptr;

		for (auto p = begin; p != end; ++p, ++q)
			*q = *p;
	}

	matrix_1d<T>& operator=(const matrix_1d<T>& x)
	{
		if (this != &x)
		{
			release();
			n_elems = x.n_elems;
			allocate();

			std::memcpy(ptr, x.ptr, n_elems * sizeof(T));
		}

		return *this;
	}

	matrix_1d<T>& operator=(matrix_1d<T>&& x)
	{
		if (this != &x)
		{
			release();

			n_elems = x.n_elems;
			ptr = x.ptr;

			x.n_elems = 0;
			x.ptr = nullptr;
		}

		return *this;
	}

	void clear()
	{
		release();
	}

	bool empty()
	{
		return ptr == nullptr;
	}

	size_t elems()	const
	{
		return n_elems;
	}

	size_t size()	const
	{
		return n_elems;
	}

	void resize(size_t _n_elems)
	{
		release();
		n_elems = _n_elems;
		allocate();
	}

	void resize(size_t _n_elems, T value)
	{
		release();
		n_elems = _n_elems;
		allocate();

		std::fill_n(ptr, n_elems, value);
	}

	T* data() const
	{
		return ptr;
	}

	T& operator()(size_t id_elem)
	{
		return ptr[id_elem];
	}

	T operator()(size_t id_elem) const
	{
		return ptr[id_elem];
	}

	T sum()	const
	{
		return std::accumulate(ptr, ptr + n_elems, (T)0);
	}

	void set_to_zero()
	{
		std::fill_n(ptr, n_elems, (T)0);
	}

	T min_coeff()	const
	{
		return *std::min_element(ptr, ptr + n_elems);
	}

	T max_coeff()	const
	{
		return *std::max_element(ptr, ptr + n_elems);
	}

	T norm() const
	{
		T r = (T)0;

		for (size_t i = 0; i < n_elems; ++i)
			r += ptr[i] * ptr[i];

		return ::sqrt(r);
	}

	matrix_1d<T>& sqrt()
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] = ::sqrt(ptr[i]);

		return *this;
	}

	matrix_1d<T>& pow2()
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] = ptr[i] * ptr[i];

		return *this;
	}

	matrix_1d<T>& abs()
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] = ptr[i] >= 0 ? ptr[i] : -ptr[i];

		return *this;
	}

	bool all_items_same() const
	{
		for (size_t i = 1; i < n_elems; ++i)
			if (ptr[i] != *ptr)
				return false;

		return true;
	}

	bool all_items_equal(const T x) const
	{
		for (size_t i = 0; i < n_elems; ++i)
			if (ptr[i] != x)
				return false;

		return true;
	}

	matrix_1d<T>& operator+=(const T x)
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] += x;

		return *this;
	}

	matrix_1d<T>& operator-=(const T x)
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] -= x;

		return *this;
	}

	matrix_1d<T>& operator*=(const T x)
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] *= x;

		return *this;
	}

	matrix_1d<T>& operator/=(const T x)
	{
		for (size_t i = 0; i < n_elems; ++i)
			ptr[i] /= x;

		return *this;
	}

	friend matrix_1d<T> operator-(const T x, const matrix_1d<T>& vec)
	{
		matrix_1d<T> res(vec);

		for (size_t i = 0; i < res.size(); ++i)
			res.ptr[i] = x - res.ptr[i];
		return res;
	}

	friend matrix_1d<T> operator/(const T x, const matrix_1d<T>& vec)
	{
		matrix_1d<T> res(vec);

		for (size_t i = 0; i < res.size(); ++i)
			res.ptr[i] = x / res.ptr[i];
		return res;
	}

	std::vector<bool> operator<(const T x) const
	{
		std::vector<bool> res(n_elems);

		for (size_t i = 0; i < n_elems; ++i)
			res[i] = ptr[i] < x;

		return res;
	}

	std::vector<bool> operator<=(const T x) const
	{
		std::vector<bool> res(n_elems);

		for (size_t i = 0; i < n_elems; ++i)
			res[i] = ptr[i] <= x;

		return res;
	}

	std::vector<bool> operator==(const T x) const
	{
		std::vector<bool> res(n_elems);

		for (size_t i = 0; i < n_elems; ++i)
			res[i] = ptr[i] == x;

		return res;
	}

	std::vector<bool> operator>=(const T x) const
	{
		std::vector<bool> res(n_elems);

		for (size_t i = 0; i < n_elems; ++i)
			res[i] = ptr[i] >= x;

		return res;
	}

	std::vector<bool> operator>(const T x) const
	{
		std::vector<bool> res(n_elems);

		for (size_t i = 0; i < n_elems; ++i)
			res[i] = ptr[i] > x;

		return res;
	}

	friend std::vector<bool> operator<(const T x, const matrix_1d<T>& vec)
	{
		return vec > x;
	}

	friend std::vector<bool> operator<=(const T x, const matrix_1d<T>& vec)
	{
		return vec >= x;
	}

	friend std::vector<bool> operator==(const T x, const matrix_1d<T>& vec)
	{
		return vec == x;
	}

	friend std::vector<bool> operator>=(const T x, const matrix_1d<T>& vec)
	{
		return vec <= x;
	}

	friend std::vector<bool> operator>(const T x, const matrix_1d<T>& vec)
	{
		return vec < x;
	}

	void print(const std::string &name)
	{
#ifdef REFRESH_MATRIX_ALLOW_PRINT
		std::cout << name << ":\n[";
		for (size_t i = 0; i < n_elems; ++i)
		{
			std::cout << ptr[i];
			if (i + 1 < n_elems)
				std::cout << "  ";
		}
		std::cout << "]\n";
#endif
	}

};

template<typename T>
matrix_1d<T> sign(const matrix_1d<T>& x)
{
	auto res(x);

	for (size_t i = 0; i < res.size(); ++i)
		if (res(i) < 0)
			res(i) = -1;
		else if (res(i) > 0)
			res(i) = 1;

	return res;
}

template<typename T>
matrix_1d<T> sqrt(const matrix_1d<T>& x)
{
	auto res(x);
	res.sqrt();
	return res;
}

template<typename T>
matrix_1d<T> pow2(const matrix_1d<T>& x)
{
	auto res(x);
	res.pow2();
	return res;
}

template<typename T>
matrix_1d<T> abs(const matrix_1d<T>& x)
{
	auto res(x);
	res.abs();
	return res;
}

template<typename T> 
matrix_1d<T> operator+(const matrix_1d<T>& vec, const T x)
{
	auto res(vec);
	res += x;
	return res;
}

template<typename T> 
matrix_1d<T> operator+(const T x, const matrix_1d<T>& vec)
{
	auto res(vec);
	res += x;
	return res;
}

template<typename T> 
matrix_1d<T> operator-(const matrix_1d<T>& vec, const T x)
{
	auto res(vec);
	res -= x;
	return res;
}

template<typename T> 
matrix_1d<T> operator-(const T x, const matrix_1d<T>& vec)
{
	auto res(vec);
	
	for(size_t i = 0; i < res.size(); ++i)
		res.ptr[i] = x - res.ptr[i];
	return res;
}

template<typename T> 
matrix_1d<T> operator*(const matrix_1d<T>& vec, const T x)
{
	auto res(vec);
	res *= x;
	return res;
}

template<typename T> 
matrix_1d<T> operator*(const T x, const matrix_1d<T>& vec)
{
	auto res(vec);
	res *= x;
	return res;
}

template<typename T> 
matrix_1d<T> operator/(const matrix_1d<T>& vec, const T x)
{
	auto res(vec);
	res /= x;
	return res;
}

template<typename T>
T dot_product(const matrix_1d<T>& a, const matrix_1d<T>& b)
{
	assert(a.size() == b.size());

#ifdef REFRESH_MATRIX_USE_OPENBLAS	
	return cblas_ddot(a.size(), a.data(), 1, b.data(), 1);
#else
	T res = (T)0;

	for (size_t i = 0; i < a.size(); ++i)
		res += a(i) * b(i);

	return res;
#endif
}

} // namespace refresh
