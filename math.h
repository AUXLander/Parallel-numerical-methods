#pragma once
#include <math.h>

template<class T, class InVa, typename Function>
T summary(InVa k, InVa N, Function operation)
{
	T value = 0;

	#pragma omp parallel for reduction(+:value)
	for (intptr_t i = k ; i < N; ++i)
	{
		value += operation(i);
	}

	return value;
}
