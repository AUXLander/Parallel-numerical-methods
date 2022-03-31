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
	size_t size_x;
	size_t size_y;

	matrix<T>* p_parent;

	T* p_start;

	bool lock{ true };

	matrix(size_t size_x, size_t size_y) :
		size_x(size_x), size_y(size_y), p_parent(nullptr), p_start(new T[size_x * size_y]{ 0 })
	{
		lock = true;
	}

	matrix(size_t size_x, size_t size_y, matrix<T>* p_parent, T* p_start) :
		size_x(size_x), size_y(size_y), p_parent(p_parent), p_start(p_start)
	{
		lock = false;
	}

	matrix(const matrix<T>& other) :
		size_x(other.size_x), size_y(other.size_y), p_parent(nullptr), p_start(new T[size_x * size_y]{ 0 })
	{
		lock = true;

		for (size_t i = 0; i < size_x; ++i)
		{
			for (size_t j = 0; j < size_y; ++j)
			{
				at(i, j) = other.at(i, j);
			}
		}
	}

	matrix(matrix&& other) :
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
		if (lock && p_start && p_parent == nullptr)
		{
			delete[] p_start;
		}
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

	matrix<T> from(size_t i, size_t j, size_t sz_x, size_t sz_y)
	{
		assert(p_parent == nullptr);

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

		return (m);
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


	void LU_decomposition(matrix<T>& L, matrix<T>& U)
	{
		assert(L.size_x == U.size_x);
		assert(L.size_y == U.size_y);
		assert(L.size_x == U.size_y);

		const auto& A = *this;

		size_t size = L.size_x;

		L(0, 0) = 1.0;
		U(0, 0) = A(0, 0);

		for (size_t j = 0; j < size; ++j)
		{
			for (size_t i = 0; i <= j; ++i)
			{
				const auto sum = summary<T, size_t>(0, i, [&](size_t k) { return L(i, k) * U(k, j); });

				auto sub = A(i, j) - sum;

				if (i == j && std::abs(sub) < 1e-10)
				{
					U(i, j) = 1.0;
				}
				else
				{
					U(i, j) = sub;
				}
			}

			auto _u = (T)1.0 / U(j, j);

			for (size_t i = 1; i < size; ++i)
			{
				const auto sum = summary<T, size_t>(0, j, [&](size_t k) { return L(i, k) * U(k, j); });

				L(i, j) = (A(i, j) - sum) * _u;
			}
		}
	}

	void LU_decomposition(matrix<T>& L, matrix<T>& U, size_t size)
	{
		assert(L.size_x == size_x);
		assert(L.size_y == size_y);

		assert(U.size_x == size_x);
		assert(U.size_y == size_y);

		assert(p_parent == nullptr);
		assert(L.p_parent == nullptr);
		assert(U.p_parent == nullptr);

		size_t rm_size = size_x - size; // remaining size
		size_t offset = 0;

		while (rm_size >= size)
		{
			auto A11 = from(offset, offset, size, size);
			auto A12 = from(offset, size, rm_size, size); // row
			auto A21 = from(size, offset, size, rm_size); // col

			auto A22 = from(size + offset, size + offset, rm_size, rm_size); // block

			auto L11 = L.from(offset, offset, size, size);
			auto L21 = L.from(offset + size, offset, size, rm_size); // col

			auto U11 = U.from(offset, offset, size, size);
			auto U12 = U.from(offset, offset + size, rm_size, size); // row

			A11.LU_decomposition(L11, U11);

			for (size_t shift = offset + size; shift < size_x; shift += size)
			{
				auto a = from(offset, shift, size, size);
				auto l = L.from(offset, offset, size, size);
				auto u = U.from(offset, shift, size, size);

				for (size_t j = 0; j < size; ++j)
				{
					for (size_t i = 0; i < size; ++i)
					{
						const auto sum = summary<T, size_t>(0, i, [&](size_t k) { return l(i, k) * u(k, j); });

						u(i, j) = a(i, j) - sum;
					}
				}
			}

			for (size_t shift = offset + size; shift < size_x; shift += size)
			{
				auto a = from(shift, offset, size, size);
				auto l = L.from(shift, offset, size, size);
				auto u = U.from(offset, offset, size, size);


				for (size_t i = 0; i < size; ++i)
				{
					for (size_t j = 0; j < size; ++j)
					{
						const auto sum = summary<T, size_t>(0, j, [&](size_t k) { return l(i, k) * u(k, j); });

						l(i, j) = (a(i, j) - sum) / u(j, j);
					}
				}
			}

			// A22 -= L21 * U12;

			for (size_t shift_y = offset + size; shift_y < size_y; shift_y += size)
			{
				const auto ml = L.from(shift_y, offset, size, size);

				for (size_t shift_x = offset + size; shift_x < size_x; shift_x += size)
				{
					const auto mu = U.from(offset, shift_x, size, size);
					auto ma = from(shift_y, shift_x, size, size);

					for (size_t i = 0; i < size; ++i)
					{
						for (size_t j = 0; j < size; ++j)
						{
							auto& cell = ma(i, j);

							for (size_t r = 0; r < size; ++r)
							{
								cell -= ml(i, r) * mu(r, j);
							}
						}
					}
				}
			}

			rm_size -= size;
			offset += size;
		}

		auto A11 = from(offset, offset, size, size);

		auto L11 = L.from(offset, offset, size, size);
		auto U11 = U.from(offset, offset, size, size);

		A11.LU_decomposition(L11, U11);

		for (int j = 0; j < size; ++j)
		{
			for (int i = 0; i < j; ++i)
			{
				L(i, j) = 0.0;
			}

			for (int i = j + 1; i < size; ++i)
			{
				U(i, j) = 0.0;
			}
		}
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

using matrix_t = matrix<double>;

void LU_Decomposition(double* A, double* L, double* U, int N)
{
	matrix_t bmA(N, N, nullptr, A);
	matrix_t bmL(N, N, nullptr, L);
	matrix_t bmU(N, N, nullptr, U);

	bmA.lock = false;
	bmL.lock = false;
	bmU.lock = false;

	for (size_t i = 0; i < N; ++i)
	{
		if (std::abs(bmA(i, i)) < 1e-10)
		{
			for (size_t j = 0; j < N; ++j)
			{
				std::swap(bmA(i, j), bmA((i + 1) % N, j));
			}
		}
	}

	matrix_t A_copy(bmA);

	size_t blsize = std::min<size_t>(120U, N);

	while (N % blsize)
	{
		--blsize;
	}

	bmA.LU_decomposition(bmL, bmU, blsize);

	auto LU = bmL * bmU;

	if (LU == A_copy)
	{
		std::cout << "\n\nprint: LU == A its ok\n";
	}
	else
	{
		// 1648742200 2231383348
		std::cout << "\n\nprint: LU != A !!!!! NOT OK\n";

		A_copy.print(std::cout);

		std::cout << "\n\nprint: LU\n";

		LU.print(std::cout);

		std::cout << "\n\nprint: L\n";

		bmL.print(std::cout);

		std::cout << "\n\nprint: U\n";

		bmU.print(std::cout);
	}
}