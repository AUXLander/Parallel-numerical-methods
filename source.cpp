#include <random>
#include <chrono>

#include "matrix.h"

void LU_decomposition(double* A, double* L, double* U, int n)
{
	matrix a(A, n, n);

	a.LU_decomposition(L, U, n);
}

void LU_decomposition(matrix_adaptor<double>& A, matrix_adaptor<double>& L, matrix_adaptor<double>& U, int n)
{
	LU_decomposition(A.data(), L.data(), U.data(), n);
}

int main()
{
	constexpr static size_t N = 3000;

	matrix A(N, N);
	matrix L(N, N);
	matrix U(N, N);

	srand(0); // no random

	for (size_t index = 0; index < N*N; ++index)
	{
		A.linear_access(index) = (double)(rand() % 100) + 1.0;
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

	auto t_norm = (L * U - A).norm();
	auto a_norm = A.norm();

	std::cout << t_norm / a_norm;


	return 0;
}