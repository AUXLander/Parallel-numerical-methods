#include <random>
#include <chrono>

#include "matrix.h"
#include "math.h"

void LU_decomposition(double* A, double* L, double* U, int n)
{
	matrix_adaptor<double> a(A, n, n);
	matrix_adaptor<double> l(L, n, n);
	matrix_adaptor<double> u(U, n, n);

	l(0, 0) = 1.0;
	u(0, 0) = a(0, 0);

	for (int i, j = 0; j < n; ++j)
	{
		#pragma omp parallel
		{
			#pragma omp for
			for (i = 0; i <= j; ++i)
			{
				const auto sum = summary<double, int>(0, i, [&](int k) { return l(i, k) * u(k, j); });

				u(i, j) = a(i, j) - sum;
			}

			#pragma omp for
			for (i = 1; i < n; ++i)
			{
				const auto sum = summary<double, int>(0, j, [&](int k) { return l(i, k) * u(k, j); });

				l(i, j) = (a(i, j) - sum) / u(j, j);
			}
		}
	}
}

void LU_decomposition(matrix_adaptor<double>& A, matrix_adaptor<double>& L, matrix_adaptor<double>& U, int n)
{
	LU_decomposition(A.matrix.get(), L.matrix.get(), U.matrix.get(), n);
}

int main()
{
	constexpr static size_t N = 1000;

	matrix_adaptor<double> A(N, N);
	matrix_adaptor<double> L(N, N);
	matrix_adaptor<double> U(N, N);

	srand(0); // no random

	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < N; ++j)
		{
			A[i][j] = (double)(rand() % 100) + 1.0;
		}
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	LU_decomposition(A, L, U, N);
	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

	std::cout << fp_ms.count() << " ms\n";

	//A.print(std::cout);

	std::cout << '\n' << '\n';

	//L.print(std::cout);

	std::cout << '\n' << '\n';

	//U.print(std::cout);

	return 0;
}