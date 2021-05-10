#pragma once

#include "cph.h"
#include <random>

class RndGenIntUni {
public:
	RndGenIntUni(const Int n, const Int m) \
		: gen(500), distr(n, m) { };
	RndGenIntUni(const RndGenIntUni &rhs) = delete;
	RndGenIntUni(RndGenIntUni &&rhs) = delete;
	RndGenIntUni & operator=(const RndGenIntUni &rhs) = delete;
	RndGenIntUni & operator=(RndGenIntUni &&rhs) = delete;
	inline void fillArr(Int *arr, const Int n);
	inline Int generate();
	~RndGenIntUni() = default;
protected:
	std::mt19937 gen;
	std::uniform_int_distribution<Int> distr;
};

inline void RndGenIntUni::fillArr(Int *arr, const Int n) {
	for (Int ind = 0; ind < n; ind++)
		arr[ind] = distr(gen);
}

inline Int RndGenIntUni::generate() {
	return distr(gen);
}
