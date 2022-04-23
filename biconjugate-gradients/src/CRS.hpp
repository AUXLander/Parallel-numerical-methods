#include <vector>

struct CRSMatrix
{
	int n; // Число строк в матрице
	int m; // Число столбцов в матрице
	int nz; // Число ненулевых элементов в разреженной матрице
	std::vector<double> val; // Массив значений матрицы по строкам
	std::vector<int> colIndex; // Массив номеров столбцов
	std::vector<int> rowPtr; // Массив индексов начала строк
};