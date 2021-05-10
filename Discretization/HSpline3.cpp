#include "HSpline3.h"

template class HSpline3<double>;
template class HSpline3<Complex>;

template <class TBas>
HSpline3<TBas>::HSpline3(const Grid &gr) : AHermitSpline<TBas>(gr) {
	deg = 3;
	basPerNode = 2;
	nFree = nBas = gr.getNPoints() * basPerNode;
	this->initRefPoints();
}

template <class TBas>
HSpline3<TBas>::HSpline3(const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight) \
		: AHermitSpline<TBas>(gr, bcTypeLeft, bcTypeRight, bcLeft, bcRight) {
	deg = 3; basPerNode = 2;
	nFree = nBas = gr.getNPoints() * basPerNode;
	if (bcTypeLeft != none) nFree--;
	if (bcTypeRight != none) nFree--;
	this->initRefPoints();
}

/*
HSpline3::HSpline3(const HSpline3 &rhs) \
	: AHermitSpline(rhs) { }

HSpline3 & HSpline3::operator=(const HSpline3 &rhs) {
	AHermitSpline::operator=(rhs);
	return *this;
}

HSpline3 & HSpline3::operator=(const ABasis &rhs) {
	const HSpline3 & rhs_ = dynamic_cast<const HSpline3 &>(rhs);
	this->operator=(rhs_);
	return *this;
}*/

template <class TBas>
TBas HSpline3<TBas>::f(const double x, const Int n) const {
	Int nInt = this->toInternalN(n);
	Int splNum = nInt % basPerNode;
	switch (splNum) {
		case 0: return phi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 1: return psi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case -1: return phi(x, lPoi[0], mPoi[0], rPoi[0]) + \
			bcLeft * psi(x, lPoi[0], mPoi[0], rPoi[0]);
		case -2: return phi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]) + \
			bcRight * psi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	}
	return 0.0;
}

template <class TBas>
TBas HSpline3<TBas>::d(const double x, const Int n) const {
	Int nInt = this->toInternalN(n);
	Int splNum = nInt % basPerNode;
	switch (splNum) {
		case 0: return dphi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 1: return dpsi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case -1: return dphi(x, lPoi[0], mPoi[0], rPoi[0]) + \
			bcLeft * dpsi(x, lPoi[0], mPoi[0], rPoi[0]);
		case -2: return dphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]) + \
			bcRight * dpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	}
	return 0.0;
}

template <class TBas>
TBas HSpline3<TBas>::dd(const double x, const Int n) const {
	Int nInt = this->toInternalN(n);
	Int splNum = nInt % basPerNode;
	switch (splNum) {
		case 0: return ddphi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case 1: return ddpsi(x, lPoi[nInt], mPoi[nInt], rPoi[nInt]);
		case -1: return ddphi(x, lPoi[0], mPoi[0], rPoi[0]) + \
			bcLeft * ddpsi(x, lPoi[0], mPoi[0], rPoi[0]);
		case -2: return ddphi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]) + \
			bcRight * ddpsi(x, lPoi[nBas - 1], mPoi[nBas - 1], rPoi[nBas - 1]);
	}
	return 0.0;
}

template <class TBas>
TBas HSpline3<TBas>::f(const double x, const AFunction<double> &func) const {
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
Complex HSpline3<TBas>::f(const double x, const AFunction<Complex> &func) const {
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
TBas HSpline3<TBas>::d(const double x, const AFunction<double> &func) const {
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
Complex HSpline3<TBas>::d(const double x, const AFunction<Complex> &func) const {
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
TBas HSpline3<TBas>::dd(const double x, const AFunction<double> &func) const {
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
Complex HSpline3<TBas>::dd(const double x, const AFunction<Complex> &func) const {
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