#pragma once
#include <assert.h>
#include <array>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>

template<class T>
struct row_adaptor
{
	const size_t size_x;
	const size_t size_y;

	T* const memory;

	size_t index_x;

	template<class...DimSizes>
	row_adaptor(T* ptr, size_t index_x, size_t size_x, size_t size_y)
		: memory(ptr), index_x(index_x), size_x(size_x), size_y(size_y)
	{;}
	
	T& operator[](size_t index_y)
	{
		assert(index_y < size_y);
		return memory[index_y * size_x + index_x];
	}

	const T& operator[](size_t index_y) const
	{
		assert(index_y < size_y);
		return memory[index_y * size_x + index_x];
	}
};


template<class T>
struct matrix_adaptor
{
	const size_t size_x;
	const size_t size_y;

	const size_t length;

	T* const memory;

	T* allocate()
	{
		return new T[length]{0};
	}
	
	template<class...DimSizes>
	matrix_adaptor(T *ptr, size_t size_x, size_t size_y)
		: size_x(size_x), size_y(size_y), length(size_x* size_y), memory(ptr)
	{;}

	template<class...DimSizes>
	matrix_adaptor(size_t size_x, size_t size_y)
		: size_x(size_x), size_y(size_y), length(size_x* size_y), memory(allocate())
	{;}

	row_adaptor<T> operator[](size_t index_x)
	{
		assert(index_x < size_x);
		return row_adaptor<T>(memory, index_x, size_x, size_y);
	}

	row_adaptor<T> operator[](size_t index_x) const
	{
		assert(index_x < size_x);
		return row_adaptor<T>(memory, index_x, size_x, size_y);
	}

	T& operator()(size_t index_y, size_t index_x)
	{
		assert(index_x < size_x);
		assert(index_y < size_y);

		return memory[index_y * size_x + index_x];
	}

	void print(std::ostream& fd)
	{
		auto stored_flags = fd.flags();

		fd << std::fixed << std::setprecision(3);

		for (size_t index = 0; index < length; ++index)
		{
			fd << memory[index] << (index % size_x == size_x - 1 ? '\n' : ' ');
		}

		fd.setf(stored_flags);
	}
};