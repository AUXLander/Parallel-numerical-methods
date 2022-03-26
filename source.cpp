#include <random>
#include <chrono>

#include "matrix.h"
#include "block_matrix.h"

void LU_decomposition(double* A, double* L, double* U, int n)
{
	matrix a(A, n, n);

	a.LU_decomposition(L, U, n);
}

int main()
{
	constexpr static size_t N = 9;
	constexpr static size_t block_size = 3;
	constexpr static size_t blocks_count = N / block_size;

	srand(0); // no random

	block_matrix<double> A(block_size, block_size, N, N);
	block_matrix<double> L(block_size, block_size, N, N);
	block_matrix<double> U(block_size, block_size, N, N);

	#pragma omp parallel 
	{
		for (int i = 0; i < blocks_count; ++i)
		{
			#pragma omp for 
			for (int j = 0; j < blocks_count; ++j)
			{
				A(i, j).fill_random();
			}
		}
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	#pragma omp parallel 
	{
		for (int i = 0; i < blocks_count; ++i)
		{
			#pragma omp for 
			for (int j = 0; j < blocks_count; ++j)
			{
				A(i, j).LU_decomposition(L(i, j).data(), U(i, j).data(), block_size);
			}
		}
	}

	auto t2 = std::chrono::high_resolution_clock::now();


	// printing
	if constexpr (true)
	{
		std::cout << "print A: \n";

		A.print(std::cout);

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

	////A.print(std::cout);

	//std::cout << '\n' << '\n';

	////L.print(std::cout);

	//std::cout << '\n' << '\n';

	////U.print(std::cout);
	//
	//return 0;

	//auto t_norm = (L * U - A).norm();
	//auto a_norm = A.norm();

	//std::cout << t_norm / a_norm;


	return 0;
}