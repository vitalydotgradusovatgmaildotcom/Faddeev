
#pragma once

#include "AnOperator.h"
#define REORTHOG_THRESH 0.25

template <class T> class IRAMEigenSolver;

template<class T>
class KrylovSubspace {
public:
	KrylovSubspace(AnOperator<T> &op, \
		const AFunction<T> &v0, const Int &nmax);
	KrylovSubspace(const KrylovSubspace<T> &rhs)  = delete;
	KrylovSubspace(KrylovSubspace<T> &&rhs) = delete;
	KrylovSubspace & operator=(const KrylovSubspace<T> &rhs) = delete;
	KrylovSubspace & operator=(KrylovSubspace<T> &&rhs) = delete;
	void makeGramSchmidtStep();
	Int dim();
	const AFunction<T> & operator[](const Int i) const;
	const GenMatrix<T> getHbar() const;
	void renew(const AFunction<T> &v0);
	~KrylovSubspace(void);
protected:
	Int n, nmax;
	AnOperator<T> &op;
	GenMatrix<T> Hbar;
	std::vector<unique_ptr<AFunction<T>>> V;
	void reduce(const Int newdim);
	friend class IRAMEigenSolver<T>;
};

