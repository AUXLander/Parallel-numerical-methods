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
		blocks(allocate(length / block_length), [](matrix* p) { delete reinterpret_cast<uint8_t*>(p); })
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