#include <random>
#include <chrono>
#include <time.h>
#include <thread>

#include "block_matrix_2.h"


int main()
{
	while (true)
	{
		unsigned int seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

		//unsigned int seed = time(NULL);

		//seed = 1648742200;

		//srand(seed = 2231383348);
		srand(seed);

		constexpr size_t N = 20;

		matrix_t A(N, N);
		matrix_t U(N, N);
		matrix_t L(N, N);

		size_t idx = 0;

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				A(i, j) = 100.0 * ((double)rand() / (double)RAND_MAX);

				++idx;
			}
		}

		

		if constexpr (N <= 10)
		{
			A.print(std::cout, false);
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
		}


	}

	return 0;
}