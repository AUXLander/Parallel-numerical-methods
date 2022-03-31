#pragma once
#include "matrix.h"

template<class T>
struct block_matrix
{
	using matrix = matrix_adaptor<T>;

	const size_t block_size_x;
	const size_t block_size_y;
	const size_t block_length;

	const size_t size_x;
	const size_t size_y;
	const size_t length;

	size_t block_count_x;
	size_t block_count_y;

	std::unique_ptr<T[], std::function<void(T*)>> memory;
	std::unique_ptr<matrix[], std::function<void(matrix*)>> blocks;

	inline matrix* allocate()
	{
		size_t x = size_x / block_size_x;
		size_t frac_x = (size_x % block_size_x > 0) ? 1 : 0;

		block_count_x = x + frac_x;

		size_t y = size_y / block_size_y;
		size_t frac_y = (size_y % block_size_y > 0) ? 1 : 0;

		block_count_y = y + frac_y;

		return reinterpret_cast<matrix*>(new uint8_t[sizeof(matrix) * (block_count_x + block_count_y)]);
	}

	block_matrix(size_t block_size_x, size_t block_size_y, size_t size_x, size_t size_y) :
		block_size_x(block_size_x), block_size_y(block_size_y), block_length(block_size_x * block_size_y),
		size_x(size_x), size_y(size_y), length(size_x * size_y),
		memory(new T[length]{ 0 }, [](T* p) { delete[] p; }),
		blocks(allocate(), [](matrix* p) { delete[] reinterpret_cast<uint8_t*>(p); })
	{
		auto* p = blocks.get();
		auto* l = memory.get();

		size_t offset = 0;
		size_t index = 0;

		for (size_t y = 0; y < size_y;)
		{
			size_t __size_y = std::min(block_size_y, size_y - y);

			for (size_t x = 0; x < size_x;)
			{
				size_t __size_x = std::min(block_size_x, size_x - x);

				new (p + index) matrix(l + __size_x * __size_y, __size_x, __size_y);

				offset += __size_x * __size_y;

				++index;

				x += __size_x;
			}

			y += __size_y;
		}
	}

	block_matrix(T* ptr, size_t block_size_x, size_t block_size_y, size_t size_x, size_t size_y) :
		block_size_x(block_size_x), block_size_y(block_size_y), block_length(block_size_x* block_size_y),
		size_x(size_x), size_y(size_y), length(size_x* size_y),
		memory(ptr, [](T* p) {;}),
		blocks(allocate(), [](matrix* p) { delete[] reinterpret_cast<uint8_t*>(p); })
	{
		auto* p = blocks.get();
		auto* l = memory.get();

		size_t offset = 0;
		size_t index = 0;

		for (size_t y = 0; y < size_y;)
		{
			size_t __size_y = std::min(block_size_y, size_y - y);

			for (size_t x = 0; x < size_x;)
			{
				size_t __size_x = std::min(block_size_x, size_x - x);

				new (p + index) matrix(l + __size_x * __size_y, __size_x, __size_y);

				offset += __size_x * __size_y;

				++index;

				x += __size_x;
			}

			y += __size_y;
		}
	}

	inline matrix& operator()(size_t i, size_t j)
	{
		const auto count_x = size_x / block_size_x;
		const auto count_y = size_y / block_size_y;

		return (blocks.get())[j + count_y * i];
	}

	operator T* ()
	{
		return memory.get();
	}

	void LU_decomposition(double* pL, double* pU, size_t block_size, size_t N)
	{
		//const size_t blocks_count = N / block_size;

		block_matrix<double>& A = *this;
		block_matrix<double> L(pL, block_size, block_size, N, N);
		block_matrix<double> U(pU, block_size, block_size, N, N);

		// k is diagonal index
		for (int k = 0; k < block_count_x; ++k)
		{
			A(k, k).LU_decomposition(L(k, k), U(k, k));
			
			for (intptr_t l = k; l < block_count_x - 1; ++l)
			{
				A(k, l + 1).krum(L(k, k), U(k, l + 1), block_size);
			}
			
			#pragma omp parallel for
			for (intptr_t l = k; l < block_count_y - 1; ++l)
			{
				A(l + 1, k).kdlm(L(l + 1, k), U(k, k), block_size);
			}

			if (k < block_count_x - 1)
			{				
				for (size_t y = k + 1; y < block_count_y; ++y)
				{
					const auto& ml = L(y, k);

					for (size_t x = k + 1; x < block_count_x; ++x)
					{
						const auto& mu = U(k, x);
						auto& ma = A(y, x);
						
					    #pragma omp for 
						for (intptr_t i = 0; i < ma.size_y; ++i)
						{
							for (size_t j = 0; j < std::min(ml.size_y, mu.size_x); ++j)
							{
								T lu = 0;

								//#pragma omp parallel for reduction(+:lu)
								for (size_t r = 0; r < std::min(ml.size_x, mu.size_y); ++r)
								{
									lu += ml(i, r) * mu(r, j);
								}

								ma(i, j) -= lu;
							}
						}
					}
				}
			}
		}
	}

	void print(std::ostream& fd)
	{
		auto stored_flags = fd.flags();

		const auto count_x = size_x / block_size_x;
		const auto count_y = size_y / block_size_y;

		fd << std::fixed << std::setprecision(3);

		for (size_t i = 0; i < count_x; ++i)
		{
			if (i > 0)
			{
				for (size_t k = 0; k < count_y * block_size_x * 8 + count_y * block_size_x * 2 - 2 + (count_y - 1) * 3; ++k)
				{
					//fd << '-';
				}

				fd << '\n';
			}

			for (size_t l = 0; l < block_size_y; ++l)
			{
				for (size_t j = 0; j < count_y; ++j)
				{
					if (j > 0)
					{
						//fd << "|  ";
						fd << "  ";
					}

					(blocks.get())[j + count_y * i].print(fd, l);
				}

				fd << '\n';
			}
		}

		fd.setf(stored_flags);
	}
};