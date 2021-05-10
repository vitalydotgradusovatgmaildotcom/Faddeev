#pragma once

#include "ABasis.h"

template <class T, Int d>
class ADiscretization {
public:
	ADiscretization();
	ADiscretization(const ADiscretization<T, d> &rhs) = delete;
	ADiscretization(ADiscretization &&rhs) = delete;
	ADiscretization<T, d> & operator=(const ADiscretization<T, d> &rhs) = delete;
	ADiscretization<T, d> & operator=(ADiscretization<T, d> &&rhs) = delete;
	virtual T f(const std::array<double, d> &xs, const AFunction<T> &func) const = 0;
	virtual T df(const std::array<double, d> &xs, const AFunction<T> &func) const = 0;
	virtual T d2f(const std::array<double, d> &xs, const AFunction<T> &func) const = 0;
	inline Int getRaw(const std::array<Int, d> &is) const;
	inline void getInds(const Int i, std::array<Int, d> &is) const;
	inline Int getN() const;
	Int getNi(const Int i) const;
	virtual ~ADiscretization(void);
protected:
	std::array<Int, d> ns = {0};
	Int n = 0;
};

template <class T, Int d>
inline Int ADiscretization<T, d>::getN() const {
	return n;
}

template <class T, Int d>
inline Int ADiscretization<T, d>::getRaw( \
		const std::array<Int, d> &is) const {
	assert(is[0] < ns[0]);
	Int res = is[0];
	for (Int k = 1; k < d; k++) {
		assert(is[k] < ns[k]);
		res = res*ns[k] + is[k];
	}
	return res;
}

template <class T, Int d>
inline void ADiscretization<T, d>::getInds( \
		const Int i, std::array<Int, d> &is) const {
	assert(i < n);
	static auto qr = std::div((Int)i, (Int)1);
	qr.quot = i;
	for (Int k = d-1; k >= 0; k--) {
		qr = std::div(qr.quot, ns[k]);
		is[k] = qr.rem;
	}
}

template <class T, Int d>
ADiscretization<T, d>::ADiscretization() { }

template <class T, Int d>
Int ADiscretization<T, d>::getNi(const Int i) const {
	return ns[i];
}

template <class T, Int d>
ADiscretization<T, d>::~ADiscretization(void) { }
