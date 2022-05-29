struct heat_task 
{
	double X; // ширина пластины
	double Y; // длина пластины
	int n; // размер сетки по x
	int m; // размер сетки по y

	double left_condition(double y) const;	// функция, задающая граничное условие при x = 0
	double right_condition(double y) const;	// функция, задающая граничное условие при x = X
	double bottom_condition(double x) const;	// функция, задающая граничное условие при y = 0
	double top_condition(double x) const;		// функция, задающая граничное условие при y = Y
	double f(double x, double y) const;		// функция, задающая внешнее воздействие
};

double heat_task::left_condition(double y) const
{
	double x = 0;

	return x*x + y * y;
}

double heat_task::right_condition(double y) const
{
	double x = 1;

	return x * x + y * y;
}

double heat_task::bottom_condition(double x) const
{
	double y = 0;

	return x * x + y * y;
}

double heat_task::top_condition(double x) const
{
	double y = 1;

	return x * x + y * y;
}

double heat_task::f(double x, double y) const
{
	return 4;
}

/////////////////////// SOFTGRADER PART BEGIN

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>

constexpr double M_PI = 3.1415;

void initialize(const heat_task& ht, 
	            double * const v, 
	            const double h,
	            const double k)
{
	const auto size_x = ht.m + 1;
	const auto size_y = ht.n + 1;

	std::fill_n(v, size_x * size_y, 0.0);

	// 

	double* vtop = v;
	double* vbottom = v + ht.m;

	#pragma omp parallel for
	for (int index = 0; index < size_y; ++index)
	{
		vtop[index * size_x] = ht.bottom_condition(h * index);
		vbottom[index * size_x] = ht.top_condition(h * index);
	}

	//

	double* vleft = v;
	double* vright = v + ht.n * size_x;

	#pragma omp parallel for
	for (int index = 0; index < size_x; ++index)
	{
		vleft[index] = ht.left_condition(k * index);
		vright[index] = ht.right_condition(k * index);
	}
}

void heat_dirichlet_sor(heat_task ht, double* v) 
{
	const auto size_x = ht.m + 1;
	const auto size_y = ht.n + 1;

	const auto h1 = ht.X / static_cast<double>(ht.n);
	const auto k1 = ht.Y / static_cast<double>(ht.m);

	const auto h2 = std::pow(h1, -2.0); 
	const auto k2 = std::pow(k1, -2.0);

	const auto omega = 2.0 / (1.0 + std::sin(M_PI * std::min(h1, k1) / 2.0));

	const auto denominator = 2.0 * (h2 + k2);

	std::unique_ptr<double[]> f{ new double[size_x * size_y] {0.0} };

	#pragma omp parallel for
	for (int i = 0; i < size_y; ++i)
	{
		for (int j = 0; j < size_x; ++j)
		{
			f[i * size_x + j] = -ht.f(i * h1, j * k1);
		}
	}

	initialize(ht, v, h1, k1);

	const int max_iteration = 1.5 / std::min(h1, k1) / M_PI * std::log(1.0 / 1e-7);

	for (int iteration = 0; iteration < max_iteration; ++iteration) 
	{
		for (int k = 0; k < ht.n + ht.m - 3; ++k) 
		{
			const int begin = std::min(1 + k, ht.n - 1);
			const int end   = std::max(1, k - ht.m + 3);

			#pragma omp parallel for
			for (int i = begin; i >= end; --i)
			{
				const auto j = ht.m - (k - i + 2);

				const auto index = i * size_x + j;

				const auto top = v[index + 1];
				const auto left = v[index - size_x];
				const auto right = v[index + size_x];
				const auto bottom = v[index - 1];

				const auto eq_h2 = h2 * (left + right);
				const auto eq_k2 = k2 * (top + bottom);

				v[index] = (omega + 0.0) * (eq_h2 + eq_k2 + f[i * size_x + j]) + \
						   (1.0 - omega) * denominator * v[index];

				v[index] /= denominator;
			}
		}
	}
}

///////////////////// SOFTGRADER PART END

int main()
{
	heat_task ht;

	ht.n = 10;
	ht.m = 10;

	ht.X = 1.0;
	ht.Y = 1.0;

	const auto h1 = ht.X / static_cast<double>(ht.n);
	const auto k1 = ht.Y / static_cast<double>(ht.m);

	std::unique_ptr<double[]> mem1(new double[(ht.n + 1) * (ht.m + 1)] {0.0});
	std::unique_ptr<double[]> mem3(new double[(ht.n + 1) * (ht.m + 1)] {0.0});

	heat_dirichlet_sor(ht, mem1.get());

	std::cout << std::setprecision(3);

	for (int i = 0; i < (ht.n + 1); ++i)
	{
		for (int j = 0; j < (ht.m + 1); ++j)
		{
			mem3[i * (ht.m + 1) + j] = std::pow(h1 * i, 2) + std::pow(k1 * j, 2);
		}
	}

	for (int i = 0; i < (ht.n + 1); ++i)
	{
		for (int j = 0; j < (ht.m + 1); ++j)
		{
			std::cout << std::setw(8) << std::fixed;

			const auto l = mem1[i * (ht.m + 1) + j];
			const auto k = mem3[i * (ht.m + 1) + j];

			const auto delta = l - k;

			if (std::abs(delta) > 0.0)
			{
				std::cout << delta  << ' ';
			}
		}

		std::cout << '\n';
	}

	return 0;
}
