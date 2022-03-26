#pragma once
#include <assert.h>
#include <array>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>
#include <memory>

#include "math.h"

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
class matrix_adaptor
{
	std::unique_ptr<T[], std::function<void(T*)>> matrix;
public:
	const size_t size_x;
	const size_t size_y;

	const size_t length;

	matrix_adaptor(matrix_adaptor<T>&& other) :
		matrix(std::move(other.matrix)),
		size_x(other.size_x), size_y(other.size_y),
		length(other.length)
	{;}

	matrix_adaptor(T *ptr, size_t size_x, size_t size_y) :
		matrix(ptr, [](T*) {}),
		size_x(size_x), size_y(size_y),
		length(size_x * size_y)
	{;}
	
	matrix_adaptor(size_t size_x, size_t size_y) : 
		matrix(new T[size_x * size_y]{0}, [](T* p) { delete[] p; }),
		size_x(size_x), size_y(size_y), length(size_x* size_y)
	{;}

	T& linear_access(size_t linear_index)
	{
		auto memory = matrix.get();

		assert(linear_index < length);

		return memory[linear_index];
	}

	row_adaptor<T> operator[](size_t index_x)
	{
		auto memory = matrix.get();

		assert(index_x < size_x);
		return row_adaptor<T>(memory, index_x, size_x, size_y);
	}

	row_adaptor<T> operator[](size_t index_x) const
	{
		auto memory = matrix.get();

		assert(index_x < size_x);
		return row_adaptor<T>(memory, index_x, size_x, size_y);
	}

	inline T& operator()(size_t index_y, size_t index_x)
	{
		auto memory = matrix.get();

		assert(index_x < size_x);
		assert(index_y < size_y);

		return memory[index_y * size_x + index_x];
	}

	inline const T& operator()(size_t index_y, size_t index_x) const
	{
		auto memory = matrix.get();

		assert(index_x < size_x);
		assert(index_y < size_y);

		return memory[index_y * size_x + index_x];
	}

	matrix_adaptor<T> operator-(const matrix_adaptor<T>& other) const
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		matrix_adaptor<T> temp(size_x, size_y);

		auto a = matrix.get();
		auto b = other.matrix.get();
		auto c = temp.matrix.get();

		for (size_t index = 0; index < length; ++index)
		{
			c[index] = a[index] - b[index];

			if (c[index] != 0.0)
			{
				c[index] = c[index];
			}
		}
		
		return temp;
	}

	matrix_adaptor<T>& operator-=(const matrix_adaptor<T>& other)
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		auto a = matrix.get();
		auto b = other.matrix.get();
		auto c = matrix.get();

		for (size_t index = 0; index < length; ++index)
		{
			c[index] = a[index] - b[index];

			if (c[index] != 0.0)
			{
				c[index] = c[index];
			}
		}

		return *this;
	}

	matrix_adaptor<T> operator+(const matrix_adaptor<T>& other) const
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		matrix_adaptor<T> temp(size_x, size_y);

		auto& a = *matrix;
		auto& b = *other.matrix;
		auto& c = *temp.matrix;

		for (size_t index = 0; index < length; ++index)
		{
			c[index] = a[index] + b[index];
		}

		return temp;
	}

	matrix_adaptor<T>& operator*=(const matrix_adaptor<T>& other)
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		//matrix_adaptor<T> temp(size_x, size_y);

		auto& a = *this;
		auto& b = other;
		auto& c = *this;

		for (size_t i = 0; i < size_x; ++i)
		{
			for (size_t j = 0; j < size_y; ++j)
			{
				c(i,j) = summary<double, size_t>(0U, size_x, [&](size_t r) { return a(i,r) * b(r,j); });
			}
		}

		return *this;
	}

	matrix_adaptor<T> operator*(const matrix_adaptor<T>& other) const
	{
		//assert(size_x == other.size_x);
		//assert(size_y == other.size_y);

		matrix_adaptor<T> temp(other.size_x, size_y);

		auto& a = *this;
		auto& b = other;
		auto& c = temp;

		for (size_t i = 0; i < other.size_x; ++i)
		{
			for (size_t j = 0; j < size_y; ++j)
			{
				c(i, j) = summary<double, size_t>(0U, size_x, [&](size_t r) { return a(i, r) * b(r, j); });
			}
		}

		return temp;
	}

	void krum(int j, double* L, double* U, int n)
	{
		const auto& a = *this;

		matrix_adaptor<T> l(L, n, n);
		matrix_adaptor<T> u(U, n, n);

		for (int z = 0; z < n; ++z)
			for (int i = 0; i < n; ++i)
			{
				const auto sum = summary<double, int>(0, i, [&](int k) { return l(i, k) * u(k, z); });

				u(i, z) = a(i, z) - sum;
			}
	}

	void kdlm(int i, double* L, double* U, int n)
	{
		const auto& a = *this;

		matrix_adaptor<T> l(L, n, n);
		matrix_adaptor<T> u(U, n, n);

		for (int z = 0; z + i < n; ++z)
			for (int j = 0; j < n; ++j)
			{
				const auto sum = summary<double, int>(0, j, [&](int k) { return l(i + z, k) * u(k, j); });

				l(i + z, j) = (a(i + z, j) - sum) / u(j, j);
			}
	}


	void LU_decomposition(double* L, double* U, int n)
	{
		const auto& a = *this;

		matrix_adaptor<T> l(L, n, n);
		matrix_adaptor<T> u(U, n, n);

		l(0, 0) = 1.0;
		u(0, 0) = a(0, 0);

		for (int i, j = 0; j < n; ++j)
		{
			for (i = 0; i <= j; ++i)
			{
				const auto sum = summary<double, int>(0, i, [&](int k) { return l(i, k) * u(k, j); });

				u(i, j) = a(i, j) - sum;
			}

			for (i = 1; i < n; ++i)
			{
				const auto sum = summary<double, int>(0, j, [&](int k) { return l(i, k) * u(k, j); });

				l(i, j) = (a(i, j) - sum) / u(j, j);
			}
		}
	}

	void fill_random()
	{
		//#pragma omp parallel for
		for (intptr_t index = 0; index < length; ++index)
		{
			linear_access(index) = (double)(rand() % 100) + 1.0;
		}
	}

	void print(std::ostream& fd)
	{
		auto memory = matrix.get();
		auto stored_flags = fd.flags();

		fd << std::fixed << std::setprecision(3) ;

		for (size_t index = 0; index < length; ++index)
		{
			fd << std::setw(6) << memory[index] << (index % size_x == size_x - 1 ? "\n" : "   ");
		}

		fd.setf(stored_flags);
	}

	void print(std::ostream& fd, size_t i)
	{
		auto memory = matrix.get();
		auto stored_flags = fd.flags();
		auto offset = size_x * i;

		fd << std::fixed << std::setprecision(3);

		for (size_t index = 0; index < size_x; ++index)
		{
			fd << std::setw(8) << memory[offset + index] << "  ";
		}

		fd.setf(stored_flags);
	}

	T norm() const
	{
		auto s = summary<double, size_t>(0U, length, [p = matrix.get()](size_t i) { return p[i] * p[i]; });

		return std::sqrt(s);
	}

	T* data() const
	{
		return matrix.get();
	}
};

using matrix = matrix_adaptor<double>;