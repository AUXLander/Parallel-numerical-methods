#pragma once
#include "math.h"
#include <iostream>
#include <memory>
#include <functional>
#include <iomanip>
#include <assert.h>

template<class T>
struct matrix
{
	constexpr static bool ENABLE_PA_LU_SWAP = false;

	size_t size_x;
	size_t size_y;

private:

	matrix<T>* p_parent;

	T* p_start;

	mutable bool enable_memory_deallocation{ true };

	T* allocation(size_t size)
	{
		return new T[size]{ (T)0.0 };
	}

public:

	matrix(size_t size_x, size_t size_y) :
		size_x(size_x), size_y(size_y),
		p_parent(nullptr), p_start(allocation(size_x* size_y))
	{
		;
	}

	matrix(size_t size_x, size_t size_y, matrix<T>* p_parent, T* p_start) :
		size_x(size_x), size_y(size_y),
		p_parent(p_parent), p_start(p_start)
	{
		prevent_memory_deallocation();
	}

	matrix(const matrix<T>& other) :
		size_x(other.size_x), size_y(other.size_y),
		p_parent(nullptr), p_start(allocation(size_x* size_y))
	{
#pragma omp parallel for
		for (intptr_t i = 0; i < size_x; ++i)
		{
			for (intptr_t j = 0; j < size_y; ++j)
			{
				at(i, j) = other.at(i, j);
			}
		}
	}

	matrix(matrix&& other) noexcept :
		size_x(other.size_x), size_y(other.size_y),
		p_parent(other.p_parent), p_start(other.p_start)
	{
		other.size_x = 0;
		other.size_y = 0;
		other.p_parent = nullptr;
		other.p_start = nullptr;
	}

	~matrix()
	{
		if (enable_memory_deallocation && !is_child() && p_start)
		{
			delete[] p_start;
		}
	}

	void prevent_memory_deallocation() const
	{
		enable_memory_deallocation = false;
	}

	bool is_child() const
	{
		return p_parent != nullptr;
	}

	operator T* ()
	{
		return p_start;
	}

	T& at(size_t i, size_t j)
	{
		assert(i < size_y);
		assert(j < size_x);

		size_t size = p_parent ? p_parent->size_x : size_x;

		return *(p_start + size * i + j);
	}

	const T& at(size_t i, size_t j) const
	{
		assert(i < size_y);
		assert(j < size_x);

		size_t size = p_parent ? p_parent->size_x : size_x;

		return *(p_start + size * i + j);
	}

	T& operator()(size_t i, size_t j)
	{
		return at(i, j);
	}

	const T& operator()(size_t i, size_t j) const
	{
		return at(i, j);
	}

	matrix<T> block(size_t i, size_t j, size_t sz_x, size_t sz_y)
	{
		assert(!is_child());

		return matrix<T>(sz_x, sz_y, this, p_start + i * size_x + j);
	}

	matrix<T> operator*(const matrix<T>& other) const
	{
		assert(other.size_x == size_y);

		matrix<T> C(other.size_x, size_y);

		const auto& a = *this;
		const auto& b = other;

		for (size_t i = 0; i < size_y; ++i)
		{
			for (size_t j = 0; j < size_y; ++j)
			{
				C(i, j) = summary<T, size_t>(0U, size_x, [&](size_t r) { return a.at(i, r) * b.at(r, j); });
			}
		}

		return C;
	}

	matrix<T>& operator-=(const matrix<T>& other)
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		for (size_t i = 0; i < size_y; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				at(i, j) -= other.at(i, j);
			}
		}

		return *this;
	}

	matrix<T> operator-(const matrix<T>& other) const
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		matrix<T> m(*this);

		m -= other;

		return m;
	}


	bool operator==(const matrix<T>& other) const
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		for (size_t i = 0; i < size_y; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				if (std::abs(at(i, j) - other.at(i, j)) > 1e-10)
				{
					return false;
				}
			}
		}

		return true;
	}


	bool operator!=(const matrix<T>& other) const
	{
		assert(size_x == other.size_x);
		assert(size_y == other.size_y);

		for (size_t i = 0; i < size_y; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				if (std::abs(at(i, j) - other.at(i, j)) > 1e-15)
				{
					return true;
				}
			}
		}

		return false;
	}

	void LU_decomposition(matrix<T>& L, matrix<T>& U) const
	{
		assert(L.size_x == U.size_x);
		assert(L.size_y == U.size_y);
		assert(L.size_x == U.size_y);

		const auto& A = *this;

		size_t size = L.size_x;

		L(0, 0) = (T)1.0;
		U(0, 0) = A(0, 0);

		for (intptr_t j = 0U; j < size; ++j)
		{
			for (intptr_t i = 0; i <= j; ++i)
			{
				const auto sum = summary<T, size_t>(0, i, [&](size_t k) { return L(i, k) * U(k, j); });

				auto sub = A(i, j) - sum;

				if (i == j && std::abs(sub) < 1e-10)
				{
					U(i, j) = (T)1.0;
				}
				else
				{
					U(i, j) = sub;
				}
			}

			for (intptr_t i = 1; i < size; ++i)
			{
				const auto sum = summary<T, size_t>(0, j, [&](size_t k) { return L(i, k) * U(k, j); });

				L(i, j) = (A(i, j) - sum) / U(j, j);
			}
		}
	}

	// block LU decomposition
	void LU_decomposition(matrix<T>& L, matrix<T>& U, size_t size) const
	{
		assert(L.size_x == size_x);
		assert(L.size_y == size_y);

		assert(U.size_x == size_x);
		assert(U.size_y == size_y);

		assert(!is_child());
		assert(!L.is_child());
		assert(!U.is_child());

		matrix<T> A(*this);

		size_t rm_size = size_x - size; // remaining size
		size_t offset = 0U;

		if constexpr (ENABLE_PA_LU_SWAP)
		{
			for (size_t i = 0; i < size; ++i)
			{
				if (std::abs(A(i, i)) < 1e-10)
				{
					for (size_t j = 0; j < size; ++j)
					{
						std::swap(A(i, j), A((i + 1) % size, j));
					}
				}
			}
		}

		while (rm_size >= size)
		{
			auto A11 = A.block(offset, offset, size, size);
			auto A12 = A.block(offset, size, rm_size, size); // row
			auto A21 = A.block(size, offset, size, rm_size); // col

			auto A22 = A.block(size + offset, size + offset, rm_size, rm_size); // block

			auto L11 = L.block(offset, offset, size, size);
			auto L21 = L.block(offset + size, offset, size, rm_size); // col

			auto U11 = U.block(offset, offset, size, size);
			auto U12 = U.block(offset, offset + size, rm_size, size); // row

			A11.LU_decomposition(L11, U11);

			#pragma omp parallel
			{
				for (intptr_t shift = offset + size; shift < size_x; shift += size)
				{
					auto a = A.block(offset, shift, size, size);
					auto l = L.block(offset, offset, size, size);
					auto u = U.block(offset, shift, size, size);

					intptr_t i;
					#pragma omp for private(i)
					for (intptr_t j = 0; j < size; ++j)
					{
						for (i = 0; i < size; ++i)
						{
							const auto sum = summary<T, size_t>(0, i, [&](size_t k) { return l(i, k) * u(k, j); });

							u(i, j) = a(i, j) - sum;
						}
					}
				}

				for (intptr_t shift = offset + size; shift < size_x; shift += size)
				{
					auto a = A.block(shift, offset, size, size);
					auto l = L.block(shift, offset, size, size);
					auto u = U.block(offset, offset, size, size);

					intptr_t j;
					#pragma omp for private(j)
					for (intptr_t i = 0; i < size; ++i)
					{
						for (j = 0; j < size; ++j)
						{
							const auto sum = summary<T, size_t>(0, j, [&](size_t k) { return l(i, k) * u(k, j); });

							l(i, j) = (a(i, j) - sum) / u(j, j);
						}
					}
				}
			}

			// A22 -= L21 * U12;
			#pragma omp parallel for 
			for (intptr_t shift_y = offset + size; shift_y < size_y; shift_y += size)
			{
				const auto l = L.block(shift_y, offset, size, size);

				#pragma omp parallel for 
				for (intptr_t shift_x = offset + size; shift_x < size_x; shift_x += size)
				{
					const auto u = U.block(offset, shift_x, size, size);
					auto a = A.block(shift_y, shift_x, size, size);

					for (intptr_t i = 0; i < size; ++i)
					{
						#pragma omp parallel for 
						for (intptr_t j = 0; j < size; ++j)
						{
							auto& cell = a(i, j);

							for (intptr_t r = 0; r < size; ++r)
							{
								cell -= l(i, r) * u(r, j);
							}
						}
					}
				}
			}

			rm_size -= size;
			offset += size;
		}

		auto A11 = A.block(offset, offset, size, size);
		auto L11 = L.block(offset, offset, size, size);
		auto U11 = U.block(offset, offset, size, size);

		A11.LU_decomposition(L11, U11);

		intptr_t i;
		#pragma omp parallel for private(i)
		for (intptr_t j = 0; j < size_x; ++j)
		{
			L(j, j) = (T)1.0;

			for (i = 0; i < j; ++i)
			{
				L(i, j) = (T)0.0;
			}

			for (i = j + 1; i < size_y; ++i)
			{
				U(i, j) = (T)0.0;
			}
		}
	}

	T norm() const
	{
		T value = (T)0.0;
		for (size_t i = 0U; i < size_y; ++i)
		{
			for (size_t j = 0U; j < size_x; ++j)
			{
				value += at(i, j) * at(i, j);
			}
		}

		return std::sqrt(value);
	}

	void print(std::ostream& fd, bool skip_zeros = true)
	{
		auto stored_flags = fd.flags();

		fd << std::fixed << std::setprecision(3);

		for (size_t i = 0; i < size_y; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				if (skip_zeros && std::abs(at(i, j)) < 1e-7)
				{
					fd << std::setw(8) << "        " << "  ";
				}
				else
				{
					fd << std::setw(8) << at(i, j) << "  ";
				}
			}

			fd << '\n';
		}

		fd.setf(stored_flags);
	}
};
