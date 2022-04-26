#include "CRS.hpp"
#include <time.h>
#include <iostream>
#include <memory>
#include <chrono>

template<class T, class InVa, typename Function>
inline T summary(InVa k, InVa N, Function operation)
{
	T value = 0;

	for (intptr_t i = k; i < N; ++i)
	{
		value += operation(i);
	}

	return value;
}

template<class T>
inline T sqr(T v)
{
	return v * v;
}


void crs_multiplication(CRSMatrix& A, double* x, double* b)
{
	intptr_t j, length = A.n;

	//double sub;

	for (intptr_t i = 0; i < length; ++i)
	{
		//sub = 0.0;

		for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k)
		{
			j = A.colIndex[k];

			b[i] += A.val[k] * x[j];

			if (i != j) {
				//prod[j] += a * xOld[i];
			}
		}

		//b[i] += sub;
	}
}


void crs_transpose_multiplication(CRSMatrix& At, double* x, double* b)
{
	intptr_t j, length = At.n;

	//double sub;

	for (intptr_t i = 0; i < length; ++i)
	{
		//sub = 0.0;

		for (int k = At.rowPtr[i]; k < At.rowPtr[i + 1]; ++k)
		{
			j = At.colIndex[k];

			b[j] += At.val[k] * x[i];

			if (i != j) {
				//prod[j] += a * xOld[i];
			}
		}

		//b[i] += sub;
	}
}

void SLE_Solver_CRS_BICG(CRSMatrix& A, double* b, double eps, int max_iter, double* x, int& count)
{
	int i, j, length = A.n;
	double rOC, rNum, numerator, norm, sum, a;

	std::unique_ptr<double[]> pOld{ new double[length] {0} };
	std::unique_ptr<double[]> rOld{ new double[length] {0} };

	std::unique_ptr<double[]> sOld{ new double[length] {0} };
	std::unique_ptr<double[]> zOld{ new double[length] {0} };

	std::unique_ptr<double[]> xOld{ new double[length] {1} };

	std::unique_ptr<double[]> r{ new double[length] {0} };
	std::unique_ptr<double[]> p{ new double[length] {0} };
	std::unique_ptr<double[]> s{ new double[length] {0} };
	std::unique_ptr<double[]> z{ new double[length] {0} };

	std::unique_ptr<double[]> buffer_1{ new double[length] {0} };
	std::unique_ptr<double[]> buffer_2{ new double[length] {0} };

	crs_multiplication(A, xOld.get(), buffer_1.get());

	// r^0 = b - Ax^0
	for (i = 0; i < length; ++i)
	{
		sOld[i] = zOld[i] = pOld[i] = rOld[i] = b[i] - buffer_1[i];
	}

	double slope, alpha_denominator, beta_nominator = 0.0;
	count = 0;

	while (count < max_iter)
	{
		slope = summary<double>(0, length, [&pOld, &rOld](int k) { return pOld[k] * rOld[k]; });

		for (i = 0; i < length; ++i)
		{
			buffer_1[i] = 0;
			buffer_2[i] = 0;
		}

		// A z^(k-1)
		crs_multiplication(A, zOld.get(), buffer_1.get());

		// alpha denominator
		alpha_denominator = summary<double>(0, length, [&sOld, &buffer_1](int k) { return sOld[k] * buffer_1[k]; });

		// A^T s
		// transpose multiplication works ok
		crs_transpose_multiplication(A, sOld.get(), buffer_2.get());

		// 2, 3, 4
		for (i = 0; i < length; ++i)
		{
			x[i] = xOld[i] + (slope * zOld[i]) / alpha_denominator;
			r[i] = rOld[i] - (slope * buffer_1[i]) / alpha_denominator;
			p[i] = pOld[i] - (slope * buffer_2[i]) / alpha_denominator;
		}

		// 5

		beta_nominator = summary<double>(0, length, [&p, &r](int k) { return p[k] * r[k]; });

		for (i = 0; i < length; ++i)
		{
			z[i] = r[i] + (beta_nominator * zOld[i]) / slope;
			s[i] = p[i] + (beta_nominator * sOld[i]) / slope;
		}

		norm = std::sqrt(summary<double>(0, length, [&x, &xOld](int k) { return sqr(x[k] - xOld[k]); }));

		++count;

		if ((norm < eps) || (beta_nominator == 0))
		{
			break;
		}

		// TODO: make cache optimization by creating struct
		for (i = 0; i < length; ++i)
		{
			xOld[i] = x[i];
			rOld[i] = r[i];
			pOld[i] = p[i];
			zOld[i] = z[i];
			sOld[i] = s[i];
		}
	}

	double n = 0.0;
	double q = 0.0;

	for (i = 0; i < length; ++i)
	{
		buffer_1[i] = 0;
	}

	for (i = 0; i < length; ++i)
	{
		for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k)
		{
			j = A.colIndex[k];
			a = A.val[k];

			q += sqr(a);

			buffer_1[i] += a * xOld[j];
		}
	}


	for (i = 0; i < length; ++i)
	{
		n += sqr(buffer_1[i] - b[i]);
	}

	if (auto aws = std::sqrt(n) / std::sqrt(q) > 0.01)
	{
		std::cout << aws << '\n';
	}

	if (count == max_iter)
	{
		return;
	}

	return;
}

//constexpr static size_t N = 1083;
constexpr static size_t N = 10;

int64_t seed;

int main()
{
	matrix<double> m(N, N);


	double b[N]{ 1 };// , 2, 3
	double x[N]{ 0 };// , 0, 0 };

	while (true)
	{
		srand(seed = std::chrono::high_resolution_clock::now().time_since_epoch().count());

		//srand(1650999764);
		//srand(46395394657400);
		//srand(47868953715300);

		//srand(52972159043100);

		//srand(53376133835800);

		for (size_t i = 0; i < N; ++i)
		{
			for (size_t j = 0; j < N; ++j)
			{
				auto v = 100 * (double)rand() / (double)RAND_MAX;

				m.at(i, j) = v;
			}

			b[i] = 100 * (double)rand() / (double)RAND_MAX;
		}

		CRSMatrix crs;

		m.make_crs(crs);

		int count = 0;

		SLE_Solver_CRS_BICG(crs, b, 1e-7, 10000, x, count);
	}



	//std::cout << "val: ";
	//for (size_t i = 0; i < crs.nz; ++i)
	//{
	//	std::cout << crs.val[i] << ' ';
	//}
	//std::cout << '\n';

	//std::cout << "colIndex: ";
	//for (size_t i = 0; i < crs.nz; ++i)
	//{
	//	std::cout << crs.colIndex[i] << ' ';
	//}
	//std::cout << '\n';

	//std::cout << "rowPtr: ";
	//for (size_t i = 0; i < crs.n+1; ++i)
	//{
	//	std::cout << crs.rowPtr[i] << ' ';
	//}
	//std::cout << '\n';



	for (size_t i = 0; i < N; ++i)
	{
		std::cout << x[i] << ' ';
	}

	return 0;
}