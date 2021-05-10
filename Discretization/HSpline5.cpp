#include "HSpline5.h"

template class HSpline5<double>;
template class HSpline5<Complex>;

template <class TBas>
HSpline5<TBas>::HSpline5(const Grid &gr) : AHermitSpline<TBas>(gr) {
	deg = 5;
	basPerNode = 3;
	nFree = nBas = gr.getNPoints() * basPerNode;
	this->initRefPoints();
}

template <class TBas>
HSpline5<TBas>::HSpline5(const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight) \
		: AHermitSpline<TBas>(gr, bcTypeLeft, bcTypeRight, bcLeft, bcRight) {
	deg = 5; basPerNode = 3;
	nFree = nBas = gr.getNPoints() * basPerNode;
	if (bcTypeLeft != none) nFree--;
	if (bcTypeRight != none) nFree--;
	this->initRefPoints();
}

/*
HSpline5::HSpline5(const HSpline5 &rhs) \
	: AHermitSpline(rhs) { }

HSpline5 & HSpline5::operator=(const HSpline5 &rhs) {
	AHermitSpline::operator=(rhs);
	return *this;
}

HSpline5 & HSpline5::operator=(const ABasis &rhs) {
	const HSpline5 & rhs_ = dynamic_cast<const HSpline5 &>(rhs);
	this->operator=(rhs_);
	return *this;
}*/

template <class TBas>
TBas HSpline5<TBas>::f(const double x, const Int n) const {
	Int nInt = this->toInternalN(n);
	Int splNum = nInt % basPerNode;
	switch (splNum) {
		case 0: return phi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 1: return psi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 2: return chi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case -1: return phi(x, lPoi[0], mPoi[0], rPoi[0]) + \
			bcLeft * psi(x, lPoi[0], mPoi[0], rPoi[0]);
		case -2: return phi(x, lPoi[nBas-1], mPoi[nBas - 1], rPoi[nBas - 1]) + \
			bcRight * psi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	}
	return 0.0;
}

template <class TBas>
TBas HSpline5<TBas>::d(const double x, const Int n) const {
	Int nInt = this->toInternalN(n);
	Int splNum = nInt % basPerNode;
	switch (splNum) {
		case 0: return dphi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 1: return dpsi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 2: return dchi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case -1: return dphi(x, lPoi[0], mPoi[0], rPoi[0]) + \
			bcLeft * dpsi(x, lPoi[0], mPoi[0], rPoi[0]);
		case -2: return dphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]) + \
			bcRight * dpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	}
	return 0.0;
}

template <class TBas>
TBas HSpline5<TBas>::dd(const double x, const Int n) const {
	Int nInt = this->toInternalN(n);
	Int splNum = nInt % basPerNode;
	switch (splNum) {
		case 0: return ddphi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 1: return ddpsi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 2: return ddchi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case -1: return ddphi(x, lPoi[0], mPoi[0], rPoi[0]) + \
			bcLeft * ddpsi(x, lPoi[0], mPoi[0], rPoi[0]);
		case -2: return ddphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]) + \
			bcRight * ddpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	}
	return 0.0;
}

template <class TBas>
TBas HSpline5<TBas>::f(const double x, const AFunction<double> &func) const {
	TBas res = TBas();

	NumbersList nl;
	this->getNonzero(x, nl);
	
	NumbersList::const_iterator nli;
	for (nli = nl.begin(); \
			nli != nl.end(); nli++)
		res += func[*nli]*f(x, *nli);
	if (bcTypeLeft == Dir) res += bcLeft*phi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeLeft == Neu) res += bcLeft*psi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeRight == Dir) res += bcRight*phi(x, lPoi[nBas-1], mPoi[nBas-1], rPoi[nBas-1]);
	if (bcTypeRight == Neu) res += bcRight*psi(x, lPoi[nBas-1], mPoi[nBas-1], rPoi[nBas-1]);

	return res;
}

template <class TBas>
Complex HSpline5<TBas>::f(const double x, const AFunction<Complex> &func) const {
	Complex res = zzero;

	NumbersList nl;
	this->getNonzero(x, nl);
	
	NumbersList::const_iterator nli;
	for (nli = nl.begin(); \
			nli != nl.end(); nli++)
		res += func[*nli]*f(x, *nli);
	if (bcTypeLeft == Dir) res += bcLeft*phi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeLeft == Neu) res += bcLeft*psi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeRight == Dir) res += bcRight*phi(x, lPoi[nBas-1], mPoi[nBas-1], rPoi[nBas-1]);
	if (bcTypeRight == Neu) res += bcRight*psi(x, lPoi[nBas-1], mPoi[nBas-1], rPoi[nBas-1]);

	return res;
}

template <class TBas>
TBas HSpline5<TBas>::d(const double x, const AFunction<double> &func) const {
	TBas res = TBas();

	NumbersList nl;
	this->getNonzero(x, nl);

	NumbersList::const_iterator nli;
	for (nli = nl.begin(); \
		nli != nl.end(); nli++)
		res += func[*nli] * d(x, *nli);
	if (bcTypeLeft == Dir) res += bcLeft * dphi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeLeft == Neu) res += bcLeft * dpsi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeRight == Dir) res += bcRight * dphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	if (bcTypeRight == Neu) res += bcRight * dpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);

	return res;
}

template <class TBas>
Complex HSpline5<TBas>::d(const double x, const AFunction<Complex> &func) const {
	Complex res = zzero;

	NumbersList nl;
	this->getNonzero(x, nl);

	NumbersList::const_iterator nli;
	for (nli = nl.begin(); \
		nli != nl.end(); nli++)
		res += func[*nli] * d(x, *nli);
	if (bcTypeLeft == Dir) res += bcLeft * dphi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeLeft == Neu) res += bcLeft * dpsi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeRight == Dir) res += bcRight * dphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	if (bcTypeRight == Neu) res += bcRight * dpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);

	return res;
}

template <class TBas>
TBas HSpline5<TBas>::dd(const double x, const AFunction<double> &func) const {
	TBas res = TBas();

	NumbersList nl;
	this->getNonzero(x, nl);

	NumbersList::const_iterator nli;
	for (nli = nl.begin(); \
		nli != nl.end(); nli++)
		res += func[*nli] * dd(x, *nli);
	if (bcTypeLeft == Dir) res += bcLeft * ddphi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeLeft == Neu) res += bcLeft * ddpsi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeRight == Dir) res += bcRight * ddphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	if (bcTypeRight == Neu) res += bcRight * ddpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);

	return res;
}

template <class TBas>
Complex HSpline5<TBas>::dd(const double x, const AFunction<Complex> &func) const {
	Complex res = zzero;

	NumbersList nl;
	this->getNonzero(x, nl);

	NumbersList::const_iterator nli;
	for (nli = nl.begin(); \
		nli != nl.end(); nli++)
		res += func[*nli] * dd(x, *nli);
	if (bcTypeLeft == Dir) res += bcLeft * ddphi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeLeft == Neu) res += bcLeft * ddpsi(x, lPoi[0], mPoi[0], rPoi[0]);
	if (bcTypeRight == Dir) res += bcRight * ddphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	if (bcTypeRight == Neu) res += bcRight * ddpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);

	return res;
}