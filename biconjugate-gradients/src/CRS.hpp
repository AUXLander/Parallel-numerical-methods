#include <vector>
#include <cassert>

struct CRSMatrix
{
	int n; // Число строк в матрице
	int m; // Число столбцов в матрице
	int nz; // Число ненулевых элементов в разреженной симметричной матрице, лежащих не ниже главной диагонали
	double* val; // Массив значений матрицы по строкам
	int* colIndex; // Массив номеров столбцов
	int* rowPtr; // Массив индексов начала строк
};

template<class T>
struct matrix
{
	size_t size_x;
	size_t size_y;

private:

	T* p_start;

	T* allocation(size_t size)
	{
		return new T[size]{ (T)0.0 };
	}

public:

	matrix(size_t size_x, size_t size_y) :
		size_x(size_x), size_y(size_y),
		p_start(allocation(size_x* size_y))
	{
		;
	}

	~matrix()
	{
		if (p_start)
		{
			delete[] p_start;
		}
	}

	operator T* ()
	{
		return p_start;
	}

	T& at(size_t i, size_t j)
	{
		assert(i < size_y);
		assert(j < size_x);

		size_t size = size_x;

		return *(p_start + size * i + j);
	}

	const T& at(size_t i, size_t j) const
	{
		assert(i < size_y);
		assert(j < size_x);

		size_t size = size_x;

		return *(p_start + size * i + j);
	}

	T& operator()(size_t i, size_t j)
	{
		return at(i, j);
	}

	const T& operator()(size_t i, size_t j) const
	{
		return at(i, j);
	}

	void make_crs(CRSMatrix& m) const
	{
		std::vector<double> v;
		std::vector<int> c;

		m.rowPtr = new int[size_y + 1U]{ 0 };

		size_t count = 0;

		for (size_t i = 0; i < size_y; ++i)
		{
			for (size_t j = 0; j < size_x; ++j)
			{
				auto& e = at(i, j);
				bool is_null = std::abs(e) == static_cast<T>(0);

				if (!is_null)
				{
					c.push_back(j);
					v.push_back(e);

					++count;
				}
			}

			m.rowPtr[i + 1U] = count;
		}

		m.n = size_x;
		m.m = size_y;

		m.nz = v.size();

		m.colIndex = new int[c.size()];
		std::copy(c.data(), c.data() + c.size(), m.colIndex);

		m.val = new double[v.size()];
		std::copy(v.data(), v.data() + v.size(), m.val);
	}
};
