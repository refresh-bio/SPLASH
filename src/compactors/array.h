#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdexcept>

using namespace std;

// *****************************************************************************************
//
template <class T>
class Array
{

public:
	int getWidth() const { return width; }
	int getHeight() const { return height; }
	const std::vector<T>& getData() const { return v; }
	std::vector<T>& getData() { return v; }

	// *****************************************************************************************
	//
	Array(int width, int height) : width(width), height(height), v(width* height) {}

	// *****************************************************************************************
	//
	Array() : width(0), height(0) {}

	// *****************************************************************************************
	//
	Array(const std::vector<T>& data, int width)
		: width(width), height(data.size() / width), v(data) {

		if ((size_t)width * (size_t)height != v.size()) {
			throw std::runtime_error("Error while constructing array from vector");
		}
	}


	// *****************************************************************************************
	//
	Array(std::vector<T>&& data, int width)
		: width(width), height(data.size() / width), v(std::move(data)) {

		if ((size_t)width * (size_t)height != v.size()) {
			throw std::runtime_error("Error while constructing array from vector");
		}
	}

	// *****************************************************************************************
	//
	Array(const Array<T>& ref) : width(ref.width), height(ref.height), v(ref.v) {}

	// *****************************************************************************************
	//
	Array(Array&& rhs) : width(rhs.width), height(rhs.height), v(std::move(rhs.v)) {}

	// *****************************************************************************************
	//
	Array<T>& operator=(const Array<T>& ref) = default;

	// *****************************************************************************************
	//
	void resize(int width, int height) {
		this->width = width;
		this->height = height;
		v.resize(width * height);
	}

	// *****************************************************************************************
	//
	void add_row() {
		++height;
		v.resize(width * height);
	}

	// *****************************************************************************************
	//
	void resize(int width, int height, const T& value) {
		this->width = width;
		this->height = height;
		v.resize(width * height, value);
	}

	// *****************************************************************************************
	//
	void clear() {
		v.clear();
	}

	// *****************************************************************************************
	//
	T* operator[](const int row) { return v.data() + row * width; }

	// *****************************************************************************************
	//
	const T* operator[](const int row) const { return v.data() + row * width; }


protected:
	int width;
	int height;
	std::vector<T> v;
};
