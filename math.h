#pragma once

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
