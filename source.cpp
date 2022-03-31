#include <random>
#include <chrono>

#include "block_matrix_2.h"

void LU_Decomposition(double* A, double* L, double* U, int N)
{
	size_t block_size = 100;//  std::sqrtf(N);

	matrix<double> bmA(N, N, nullptr, A);
	matrix<double> bmL(N, N, nullptr, L);
	matrix<double> bmU(N, N, nullptr, U);

	bmA.lock = false;
	bmL.lock = false;
	bmU.lock = false;

	bmA.LU_decomposition(bmL, bmU, block_size);
}

int main()
{
	srand(0);

	constexpr size_t N = 1000;

	matrix<double> A(N, N);
	matrix<double> U(N, N);
	matrix<double> L(N, N);

	size_t idx = 0;

	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			A(i, j) = 1.0 + rand()%100;

			++idx;
		}
	}

	if constexpr (N <= 10)
	{
		A.print(std::cout);
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	LU_Decomposition(A, L, U, N);

	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

	std::cout << "\n" << fp_ms.count() << " ms\n";

	if constexpr (N <= 10)
	{
		std::cout << "\n\nprint: L\n";

		L.print(std::cout);

		std::cout << "\n\nprint: U\n";

		U.print(std::cout);

		auto LU = L * U;

		std::cout << "\n\nprint: LU\n";

		LU.print(std::cout);
	}

	return 0;
}