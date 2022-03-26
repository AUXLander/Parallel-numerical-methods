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

	std::unique_ptr<T[], std::function<void(T*)>> memory;
	std::unique_ptr<matrix[], std::function<void(matrix*)>> blocks;

	inline matrix* allocate(size_t count)
	{
		return reinterpret_cast<matrix*>(new uint8_t[sizeof(matrix) * count]);
	}

	block_matrix(size_t block_size_x, size_t block_size_y, size_t size_x, size_t size_y) :
		block_size_x(block_size_x), block_size_y(block_size_y), block_length(block_size_x * block_size_y),
		size_x(size_x), size_y(size_y), length(size_x * size_y),
		memory(new T[length]{ 0 }, [](T* p) { delete[] p; }),
		blocks(allocate(length / block_length), [](matrix* p) { delete[] reinterpret_cast<uint8_t*>(p); })
	{
		auto* p = blocks.get();
		auto* l = memory.get();

		const auto count_x = size_x / block_size_x;
		const auto count_y = size_y / block_size_y;

		for (size_t i = 0; i < count_x; ++i)
		{
			for (size_t j = 0; j < count_y; ++j)
			{
				auto offset_count = j + count_y * i;

				new (p + offset_count) matrix(l + block_length * offset_count, block_size_x, block_size_y);
			}
		}
	}

	block_matrix(T* ptr, size_t block_size_x, size_t block_size_y, size_t size_x, size_t size_y) :
		block_size_x(block_size_x), block_size_y(block_size_y), block_length(block_size_x* block_size_y),
		size_x(size_x), size_y(size_y), length(size_x* size_y),
		memory(ptr, [](T* p) {;}),
		blocks(allocate(length / block_length), [](matrix* p) { delete[] reinterpret_cast<uint8_t*>(p); })
	{
		auto* p = blocks.get();
		auto* l = memory.get();

		const auto count_x = size_x / block_size_x;
		const auto count_y = size_y / block_size_y;

		for (size_t i = 0; i < count_x; ++i)
		{
			for (size_t j = 0; j < count_y; ++j)
			{
				auto offset_count = j + count_y * i;

				new (p + offset_count) matrix(l + block_length * offset_count, block_size_x, block_size_y);
			}
		}
	}

	inline matrix& operator()(size_t i, size_t j)
	{
		const auto count_x = size_x / block_size_x;
		const auto count_y = size_y / block_size_y;

		return (blocks.get())[j + count_y * i];
	}

	inline T* data()
	{
		return memory.get();
	}

	void LU_decomposition(double* pL, double* pU, size_t block_size, size_t N)
	{
		const size_t blocks_count = N / block_size;

		block_matrix<double>& A = *this;
		block_matrix<double> L(pL, block_size, block_size, N, N);
		block_matrix<double> U(pU, block_size, block_size, N, N);

		// k is diagonal index
		for (int k = 0; k < blocks_count; ++k)
		{
			A(k, k).LU_decomposition(L(k, k).data(), U(k, k).data(), block_size);

			for (int j = k; j < blocks_count - 1; ++j)
			{
				A(k, j + 1).krum(k, L(k, k).data(), U(k, j + 1).data(), block_size);
			}

			for (int i = k; i < blocks_count - 1; ++i)
			{
				A(i + 1, k).kdlm(0, L(i + 1, k).data(), U(k, k).data(), block_size);
			}

			if (k < blocks_count - 1)
			{
				size_t t = block_size * (blocks_count - k - 1);

				matrix_adaptor<double> lr(block_size, t);
				matrix_adaptor<double> ur(t, block_size);
				matrix_adaptor<double> ar(t, t);

				for (int y = k + 1, ylr = 0; y < blocks_count; ++y, ylr += L.block_size_x)
				{
					auto& block = L(y, k);

					for (int i = 0; i < L.block_size_x; ++i)
					{
						for (int j = 0; j < L.block_size_y; ++j)
						{
							lr(i + ylr, j) = block(i, j);
						}
					}
				}

				for (int x = k + 1, xur = 0; x < blocks_count; ++x, xur += L.block_size_y)
				{
					auto& block = U(k, x);

					for (int i = 0; i < L.block_size_x; ++i)
					{
						for (int j = 0; j < L.block_size_y; ++j)
						{
							ur(i, j + xur) = block(i, j);
						}
					}
				}

				for (int x = k + 1, xar = 0; x < blocks_count; ++x, xar += A.block_size_y)
				{
					for (int y = k + 1, yar = 0; y < blocks_count; ++y, yar += A.block_size_x)
					{
						auto& block = A(y, x);

						for (int i = 0; i < A.block_size_x; ++i)
						{
							for (int j = 0; j < A.block_size_y; ++j)
							{
								ar(i + yar, j + xar) = block(i, j);
							}
						}
					}
				}

				ar -= lr * ur;

				for (int x = k + 1, xar = 0; x < blocks_count; ++x, xar += A.block_size_y)
				{
					for (int y = k + 1, yar = 0; y < blocks_count; ++y, yar += A.block_size_x)
					{
						auto& block = A(y, x);

						for (int i = 0; i < A.block_size_x; ++i)
						{
							for (int j = 0; j < A.block_size_y; ++j)
							{
								block(i, j) = ar(i + yar, j + xar);
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
					fd << '-';
				}

				fd << '\n';
			}

			for (size_t l = 0; l < block_size_y; ++l)
			{
				for (size_t j = 0; j < count_y; ++j)
				{
					if (j > 0)
					{
						fd << "|  ";
					}

					(blocks.get())[j + count_y * i].print(fd, l);
				}

				fd << '\n';
			}
		}

		fd.setf(stored_flags);
	}
};