#include <random>
#include <chrono>

#include "matrix.h"
#include "block_matrix.h"

constexpr static size_t N = 6;
constexpr static size_t block_size = 3;
constexpr static size_t blocks_count = N / block_size;

void LU_Decomposition(double* A, double* L, double* U, int n)
{
	block_matrix<double> bmA(A, block_size, block_size, N, N);
	block_matrix<double> bmL(L, block_size, block_size, N, N);
	block_matrix<double> bmU(U, block_size, block_size, N, N);

	bmA.LU_decomposition(L, U, block_size, N);
}

int main()
{
	srand(0); // no random

	block_matrix<double> A(block_size, block_size, N, N);
	block_matrix<double> L(block_size, block_size, N, N);
	block_matrix<double> U(block_size, block_size, N, N);

	//#pragma omp parallel 
	{
		for (int i = 0; i < blocks_count; ++i)
		{
			//#pragma omp for 
			for (int j = 0; j < blocks_count; ++j)
			{
				A(i, j).fill_random();
			}
		}
	}

	if constexpr (N <= 10)
	{
		std::cout << "print A: \n";
		A.print(std::cout);
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	LU_Decomposition(A, L, U, N);

	auto t2 = std::chrono::high_resolution_clock::now();

	if constexpr ( N <= 10 )
	{
		std::cout << '\n' << '\n';

		std::cout << "print L: \n";
		L.print(std::cout);

		std::cout << '\n' << '\n';

		std::cout << "print U: \n";
		U.print(std::cout);
	}

	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

	std::cout << '\n' << '\n';

	std::cout << fp_ms.count() << " ms\n";

	std::cout << '\n' << '\n';

	return 0;
}