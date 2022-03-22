#include "matrix.h"

template<class T, class InVa>
T summary(InVa k, InVa N, std::function<T(InVa)>&& operation)
{
	T value = 0;

	for (; k < N; ++k)
	{
		value += operation(k);
	}

	return value;
}

void LU_dec(double* A, double* L, double* U, int n)
{
	matrix_adaptor<double> a(A, n, n);
	matrix_adaptor<double> l(L, n, n);
	matrix_adaptor<double> u(U, n, n);

	l(0, 0) = 1.0;
	u(0, 0) = a(0, 0);

	double sum;

	for (int i, j = 0; j < n; ++j)
	{
		for (i = 0; i <= j; ++i)
		{
			sum = summary<double, int>(0, i, [&](int k) { return l(i, k) * u(k, j); });

			u(i, j) = a(i, j) - sum;
		}

		for (i = 1; i < n; ++i)
		{
			sum = summary<double, int>(0, j, [&](int k) { return l(i, k) * u(k, j); });

			l(i, j) = (a(i, j) - sum) / u(j, j);
		}
	}
}


int main()
{
	matrix_adaptor<double> A(2, 2);
	matrix_adaptor<double> L(2, 2);
	matrix_adaptor<double> U(2, 2);

	A[0][0] = 1.0;
	A[1][0] = 2.0;
	A[0][1] = 3.0;
	A[1][1] = 4.0;

	LU_dec(A.memory, L.memory, U.memory, 2);

	A.print(std::cout);

	std::cout << '\n' << '\n';

	L.print(std::cout);

	std::cout << '\n' << '\n';

	U.print(std::cout);

	return 0;
}