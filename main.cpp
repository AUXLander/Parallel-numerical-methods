#include <random>
#include <chrono>
#include <thread>

#include "matrix.h"

using matrix_t = matrix<double>;

void LU_Decomposition(double* pA, double* pL, double* pU, int N)
{
	const matrix_t A(N, N, nullptr, pA);

	matrix_t L(N, N, nullptr, pL);
	matrix_t U(N, N, nullptr, pU);

	A.prevent_memory_deallocation();
	L.prevent_memory_deallocation();
	U.prevent_memory_deallocation();

	size_t blsize = std::min<size_t>(120U, N);

	while (N % blsize)
	{
		--blsize;
	}

	if (blsize == 1U)
	{
		blsize = N;
	}

	A.LU_decomposition(L, U, blsize);
}


constexpr bool ENABLE_LOOP_TESTING = false;
constexpr bool ENABLE_RESULT_CHECK = false;

constexpr bool SKIP_ZEROS_OUTPUT = true;

constexpr size_t BLOCK_MATRIX_SIZE = 10U;

int main()
{
	do
	{
		unsigned int seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

		srand(seed);

		matrix_t A(BLOCK_MATRIX_SIZE, BLOCK_MATRIX_SIZE);
		matrix_t U(BLOCK_MATRIX_SIZE, BLOCK_MATRIX_SIZE);
		matrix_t L(BLOCK_MATRIX_SIZE, BLOCK_MATRIX_SIZE);

		

		for (size_t i = 0; i < BLOCK_MATRIX_SIZE; ++i)
		{
			for (size_t j = 0; j < BLOCK_MATRIX_SIZE; ++j)
			{
				A(i, j) = 100.0 * ((double)rand() / (double)RAND_MAX);
			}
		}


		// 

		const auto t1 = std::chrono::high_resolution_clock::now();

		LU_Decomposition(A, L, U, BLOCK_MATRIX_SIZE);

		const auto t2 = std::chrono::high_resolution_clock::now();

		// 


		std::cout << "LU Decomposition source matrix A [ " << BLOCK_MATRIX_SIZE;
		std::cout << " x " << BLOCK_MATRIX_SIZE;
		std::cout << " ] (seed = " << seed << ") has been finished in ";
		std::cout << std::chrono::duration<double, std::milli>(t2 - t1).count() << " ms\n\n";

		if constexpr (BLOCK_MATRIX_SIZE <= 10)
		{
			std::cout << "print A:\n";
			A.print(std::cout, false);

			std::cout << "\n\nprint L:\n";

			L.print(std::cout, SKIP_ZEROS_OUTPUT);

			std::cout << "\n\nprint U:\n";

			U.print(std::cout, SKIP_ZEROS_OUTPUT);
		}

		if constexpr (ENABLE_RESULT_CHECK)
		{
			matrix_t LU = L * U;

			// some elements have different
			if (LU != A)
			{
				const auto diff = (LU - A).norm() / A.norm();

				if (diff >= 0.01)
				{
					std::cout << "\n\nprint: ||LU-A|| / ||A|| >= 0.01 \n";

					return 1;
				}
				else
				{
					std::cout << "\n\nprint: LU != A but ||LU-A|| / ||A|| < 0.01 - OK \n";
				}
			}
		}
	}
	while (ENABLE_LOOP_TESTING);

	return 0;
}