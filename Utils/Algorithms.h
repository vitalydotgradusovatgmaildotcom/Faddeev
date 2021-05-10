#pragma once

#include "cph.h"

//sorts first n elements of v and w simultaneously w.r.t. comparison function comp
//stable sort
template<class RandomIt1, class RandomIt2, class Compare>
inline void sortSimult(RandomIt1 vFirst, RandomIt1 vLast, \
	RandomIt2 wFirst, RandomIt2 wLast, Compare comp) {
	const Int n = std::distance(vFirst, vLast);
	assert(std::distance(wFirst, wLast) == n);
	std::vector<Int> inds; inds.resize(n);

	std::iota(inds.begin(), inds.end(), 0);
	//std::sort(inds.begin(), inds.end(), comp);
	std::stable_sort(inds.begin(), inds.end(), comp);
	Int icurr, inew;
	for (Int i = 0; i < n; i++) {
		if (inds[i] == -1) continue;
		icurr = i;
		while (inds[icurr] != i) {
			std::swap(vFirst[icurr], vFirst[inds[icurr]]);
			std::swap(wFirst[icurr], wFirst[inds[icurr]]);
			inew = inds[icurr];
			inds[icurr] = -1;
			icurr = inew;
		}
		inds[icurr] = -1;
	}
}

template <class RandomIt, class T>
inline Int countSorted(RandomIt first, RandomIt last, const T &val) {
	auto lo = lower_bound(first, last, val);
	if (lo == last || *lo != val)
		return 0;
	auto up = upper_bound(first, last, val);
	return up - lo;
}

template <class RandomIt, class T>
inline void findIntervSorted(RandomIt first, RandomIt last, const T &val, \
	Int &ind1, Int &ind2) {
	assert((val >= *first) && (val <= *(last-1)));
	if (val == *(last - 1)) {
		ind1 = last - first - 2;
		ind2 = ind1 + 1;
		return;
	}
	auto up = upper_bound(first, last, val);
	ind1 = up - first - 1;
	ind2 = ind1 + 1;
	return;
}