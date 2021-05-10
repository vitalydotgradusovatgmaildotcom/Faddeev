//var 4 --- added basis functions on last interval
#include "HybridBasis.h"
#include "AHermitSpline.h"
#include "Algorithms.h"
#include "SimplePoly.h"

template class HybridBasis<AHermitSpline<Complex>>;

template <class BasT>
HybridBasis<BasT>::HybridBasis(unique_ptr<const ABasis<Complex>> &&bas, \
	std::vector<shared_ptr<CWFCalculator>> &cwf, \
	const std::vector<double> &pn, const std::vector<double> &y0nl, \
										const std::vector<basset> &bs) : \
		ABasis(bas->getGrid()) {
//grid
	assert(grid.getNIntervals() > 1); // >= 2 intervals are assumed
//calculate number of additional functions
	nBas1 = 0;
	for (auto bs_ : bs) {
		if (bs_.deg == 0)
			nBas1++;
		else { //bs_.deg = 1
			if (bs_.pm)
				nBas1 += 3;
			else //bs_.pm = false
				nBas1 += 2;
		}
	}
//make "bas" part
	array<shared_ptr<const ABasis<Complex>>, 1> tmp = \
		array<shared_ptr<const ABasis<Complex>>, 1>({ std::move(bas) });
	discr0 = make_shared<AProjDiscr<Complex, 1, Complex>>(tmp);
	bas0 = discr0->bases[0].get();
	//throw away basis "bas" functions associated with last point
	//nFree0 = bas0->getNCoef() - bas0->nFuncOnInterval()/2;
	//nFree0 = bas0->getNCoef();
	nFree0 = bas0->getNCoef();
	if (pn.size() > 0)
		nFree0 -= 1 + nBas1;
	//nFree0 += 1;
	numBas0.reserve(nFree0);
	//Int nFree0_ = bas0->getNCoef() - bas0->nFuncOnInterval() / 2;
	//for (Int k = 0; k < nFree0_; k++)
	//	numBas0.push_back(k);
	//for (Int k = nFree0_+2; k < bas0->getNCoef(); k++)
	//	numBas0.push_back(k);
	for (Int j = 0; j < nFree0; j++)
		numBas0.push_back(j);
	numBas0_inv.resize(bas0->getNCoef());
	for (Int j = 0; j < nFree0; j++)
		numBas0_inv[numBas0[j]] = j;
	assert(bas0->getBCTypeRight() == none);
	f0 = new Function<Complex, 1>(discr0);
	for (Int j = 0; j < bas0->getNCoef(); j++)
		if (!binary_search(numBas0.begin(), numBas0.end(), j))
			f0->coef[j] = zzero;

//make additional functions part
	assert((cwf.size() == pn.size()) && (y0nl.size() == pn.size()));
	assert(bs.size() == pn.size());
	makeAddFunc(cwf, pn, y0nl, bs);

//make basis
	nBas = nFree = nFree0 + nBas1;
	if (bas0->getBCTypeLeft() != none)
		nBas++; //NB! I assume that BC fixes 1 coefficient
	bcTypeLeft = bas0->getBCTypeLeft();
	bcTypeRight = none;
	bcLeft = bas0->getBCLeft();
	bcRight = zzero;
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->f(x, numBas0[n]);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] >= iInt)
			return rnl[n1](array<double, 1>({ x }));
		else //inl[n1] < iInt
			return vfnl(n1, x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->d(x, numBas0[n]);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] >= iInt)
			return rnl[n1].df(array<double, 1>({ x }));
		else //inl[n1] < iInt
			return dfnl(n1, x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->dd(x, numBas0[n]);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] >= iInt)
			return rnl[n1].d2f(array<double, 1>({ x }));
		else //inl[n1] < iInt
			return ddfnl(n1, x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int j = 0; j < nFree0; j++)
		f0->coef[numBas0[j]] = func.coef[j];
	Complex res = bas0->f(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * vfnl(n1, x);
		it++;
	}
	while (it != inl.end()) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * rnl[n1](array<double, 1>({ x }));
		it++;
	}

	return res;
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int j = 0; j < nFree0; j++)
		f0->coef[numBas0[j]] = func.coef[j];
	Complex res = bas0->d(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * dfnl(n1, x);
		it++;
	}
	while (it != inl.end()) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * rnl[n1].df(array<double, 1>({ x }));
		it++;
	}

	return res;
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int j = 0; j < nFree0; j++)
		f0->coef[numBas0[j]] = func.coef[j];
	Complex res = bas0->dd(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * ddfnl(n1, x);
		it++;
	}
	while (it != inl.end()) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * rnl[n1].d2f(array<double, 1>({ x }));
		it++;
	}

	return res;
}

template <class BasT>
void HybridBasis<BasT>::getNonzero(const double x, NumbersList &limits) const {
	if (x < grid.getLeftmostPoint() || x > grid.getRightmostPoint()) {
		limits.resize(0);
		return;
	}
	
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0

	//"bas" part
	bas0->getNonzero(x, limits);
	//if (iInt == (grid.getNIntervals() - 1)) { //last interval
	auto it = std::remove_if(limits.begin(), limits.end(), \
		[this](Int i) {return  !binary_search(numBas0.begin(), numBas0.end(), i);});
	limits.erase(it, limits.end()); //C++20 replace with erase_if
	for (Int k = 0; k < limits.size(); k++)
		limits[k] = numBas0_inv[limits[k]];
	//}

	//additional functions part
	for (Int j = nFree0; j < nFree; j++)
		limits.push_back(j);
}

template <class BasT>
Int HybridBasis<BasT>::nFuncOnInterval() const {
	return bas0->nFuncOnInterval() + nBas1;
}

template <>
Grid HybridBasis<AHermitSpline<Complex>>::collocGrid() const {
	const AHermitSpline<Complex> *bas0_ = \
		dynamic_cast<const AHermitSpline<Complex> *>(bas0);
	if (bas0_ == nullptr)
		assert(false);

	Int basPerNode = bas0_->nFuncOnInterval() / 2;
	Int iPoi = 0;
	//Int nPoi = grid.getNIntervals()*basPerNode + nBas1;
	//if (bas0_->getBCTypeLeft() != none) nPoi--;
	//Int nAddPoi = nFree - nPoi;
	Int nAddPoi = basPerNode;
	if (inl.size() > 0)
		nAddPoi -= 1;
	nAddPoi -= (bas0_->getBCTypeLeft() != none) ? 1 : 0;
	//nAddPoi += 1;
	vector<double> poi; poi.resize(nFree);

	Int nGaussP; double a, b, half, mid;
	vector<double> gaussP, gaussW;

	//Int inl0;
	//if (inl.size() > 0)
	//	inl0 = inl.front();
	//else
	//	inl0 = 1;
	for (Int iInt = 0; iInt < grid.getNIntervals(); iInt++) {
		nGaussP = basPerNode;
		//if (iInt == inl0 && bas0_->getBCTypeLeft() != none)
		//	nGaussP--; //bcTypeLeft = Dir or Neu or Mix
		//nGaussP += countSorted(inl.begin(), inl.end(), iInt);
		if (iInt == 0)
			nGaussP += nAddPoi / 2;
		if (iInt == grid.getNIntervals() - 1)
			nGaussP += nAddPoi - nAddPoi / 2;
		//if (iInt == grid.getNIntervals() - 1)
		//	nGaussP += nAddPoi;
		gaussP.resize(nGaussP); gaussW.resize(nGaussP);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		a = grid[iInt]; b = grid[iInt + 1];
		half = 0.5*(b - a); mid = 0.5*(b + a);
		for (Int i = 0; i < gaussP.size(); i++)
			poi[iPoi++] = mid + half * gaussP[i];
	}

	return Grid(poi);
}

template <class BasT>
GenMatrix<Complex> HybridBasis<BasT>::getOverlap() const {
	assert(false);
	GenMatrix<Complex> res(0);
	return res;
}

template <class BasT>
GenMatrix<Complex> HybridBasis<BasT>::getOverlapW(\
	std::function<double(double)> rho, const Int fdeg) const {
	assert(false);
	GenMatrix<Complex> res(0);
	return res;
}

template <class BasT>
void HybridBasis<BasT>::makeAddFunc(std::vector<shared_ptr<CWFCalculator>> &cwf, \
	const std::vector<double> &pn, const std::vector<double> &y0nl, \
					const std::vector<basset> &bs) {
	std::vector<double> pn_;
	std::vector<shared_ptr<CWFCalculator>> cwf_;
	std::vector<Int> k_;
	std::vector<bool> plus_;

	//find I_{nl}
	std::vector<Int> inl_;
	Int iInt;
	inl_.reserve(y0nl.size());
	for (auto y0 : y0nl) {
		iInt = grid.getInterval(y0);
		assert(iInt != grid.getNIntervals() - 1); //y0_{nl} cannot lie in the last interval
		inl_.push_back(iInt);
	}

	inl.reserve(nBas1);
	pn_.reserve(nBas1); cwf_.reserve(nBas1);
	k_.reserve(nBas1); plus_.reserve(nBas1);
	for (Int i = 0; i < bs.size(); i++) {
		pn_.push_back(pn[i]);
		cwf_.push_back(cwf[i]);
		k_.push_back(0);
		plus_.push_back(true);
		inl.push_back(inl_[i]);
		if (bs[i].deg > 0) {
			pn_.push_back(pn[i]);
			cwf_.push_back(cwf[i]);
			k_.push_back(1);
			plus_.push_back(true);
			inl.push_back(inl_[i]);
			if (bs[i].pm) {
				pn_.push_back(pn[i]);
				cwf_.push_back(cwf[i]);
				k_.push_back(1);
				plus_.push_back(false);
				inl.push_back(inl_[i]);
			}
		}
	}

//sort functions with respect to I_{nl}
	std::vector<Int> p; p.resize(nBas1);
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(), \
			[this](Int i, Int j) {return inl[i] < inl[j]; });
	//std::vector<Int> pinv; pinv.resize(nBas1);
	//for (Int k = 0; k < nBas1; k++)
	//	pinv[p[k]] = k;

	fnl.reserve(nBas1);
	for (Int j = 0; j < nBas1; j++)
		fnl.push_back(Fnl(cwf_[p[j]], pn_[p[j]], k_[p[j]], plus_[p[j]]));
	std::sort(inl.begin(), inl.end());

//make R_{nl} functions
	makeRnl();
}

template <>
void HybridBasis<AHermitSpline<Complex>>::makeRnl() {
	const AHermitSpline<Complex> *bas0_ = \
		dynamic_cast<const AHermitSpline<Complex> *>(bas0);
	if (bas0_ == nullptr)
		assert(false);

	rnl.reserve(nBas1);
	Int nDer = bas0_->getNDeriv();
	assert(nDer <= 2); //only second derivative supported atm
	vector<Complex> val; val.resize(2*(nDer + 1));
	vector<double> poi; poi.resize(2);
	std::fill(val.begin(), val.begin()+nDer+1, zzero);
	poi[0] = grid.getLeftmostPoint();
	for (Int j = 0; j < nBas1; j++) {
		Grid gr(grid.getLeftmostPoint(), grid[inl[j]+1], 2);
		array<shared_ptr<const ABasis<Complex>>, 1> tmp;
		tmp[0] = make_shared<SimplePoly<Complex>>(gr, bas0_->getDeg());
		shared_ptr<AProjDiscr<Complex, 1, Complex>> discr = \
			make_shared<AProjDiscr<Complex, 1, Complex>>(tmp);
		rnl.push_back(Function<Complex, 1>(discr));
		poi[1] = grid[inl[j] + 1];
		for (Int iDer = 0; iDer <= nDer; iDer++) {
			switch (iDer) {
			case 0:
				val[nDer + 1] = vfnl(j, poi[1]);
				break;
			case 1:
				val[nDer + 2] = dfnl(j, poi[1]);
				break;
			case 2:
				val[nDer + 3] = ddfnl(j, poi[1]);
				break;
			}
		}
		const SimplePoly<Complex> *sp = \
			dynamic_cast<const SimplePoly<Complex> *>(discr->bases[0].get());
		sp->interpH(rnl.back(), poi, val);
	}

}

template <class BasT>
Complex HybridBasis<BasT>::vfnl(const Int j, const double y) const {
	Complex val;
	switch (fnl[j].k) {
	case 0:
		val = fnl[j].cwf->ulp(fnl[j].pn*y);
		break;
	case 1:
		if (fnl[j].plus)
			val = fnl[j].cwf->ulp(fnl[j].pn*y) / (fnl[j].pn*y);
		else
			val = fnl[j].cwf->ulm(fnl[j].pn*y) / (fnl[j].pn*y);
		break;
	default:
		assert(false); //only degree 1 correction atm
	}
	return val;
}

template <class BasT>
Complex HybridBasis<BasT>::dfnl(const Int j, const double y) const {
	Complex val;
	double pny;
	switch (fnl[j].k) {
	case 0:
		val = fnl[j].pn*fnl[j].cwf->dulp(fnl[j].pn*y);
		break;
	case 1:
		pny = fnl[j].pn*y;
		if (fnl[j].plus)
			val = fnl[j].pn*(fnl[j].cwf->dulp(pny) - fnl[j].cwf->ulp(pny) / pny) / pny;
		else
			val = fnl[j].pn*(fnl[j].cwf->dulm(pny) - fnl[j].cwf->ulm(pny) / pny) / pny;
		break;
	default:
		assert(false); //only degree 1 correction atm
	}
	return val;
}

template <class BasT>
Complex HybridBasis<BasT>::ddfnl(const Int j, const double y) const {
	Complex val;
	double pny;
	switch (fnl[j].k) {
	case 0:
		val = fnl[j].pn*fnl[j].pn*fnl[j].cwf->ddulp(fnl[j].pn*y);
		break;
	case 1:
		pny = fnl[j].pn*y;
		if (fnl[j].plus)
			val = fnl[j].pn*fnl[j].pn*( fnl[j].cwf->ddulp(pny) + \
				2.0*( -fnl[j].cwf->dulp(pny) + fnl[j].cwf->ulp(pny) /pny) / pny ) / pny;
		else
			val = fnl[j].pn*fnl[j].pn*(fnl[j].cwf->ddulm(pny) + \
				2.0*(-fnl[j].cwf->dulm(pny) + fnl[j].cwf->ulm(pny) / pny) / pny) / pny;
		break;
	default:
		assert(false); //only degree 1 correction atm
	}
	return val;
}

template <class BasT>
HybridBasis<BasT>::~HybridBasis() {
	if (f0 != nullptr)
		delete f0;
}

//var 3 --- Rnl support on [0, y2_{nl}]
/*
#include "HybridBasis.h"
#include "AHermitSpline.h"
#include "Algorithms.h"
#include "SimplePoly.h"

template class HybridBasis<AHermitSpline<Complex>>;

template <class BasT>
HybridBasis<BasT>::HybridBasis(unique_ptr<const ABasis<Complex>> &&bas, \
	std::vector<unique_ptr<CWFCalculator>> &&cwf, \
	const std::vector<double> &pn, const std::vector<double> &y0nl) : \
		ABasis(bas->getGrid()) {
//grid
	assert(grid.getNIntervals() > 1); // >= 2 intervals are assumed
//make "bas" part
	array<shared_ptr<const ABasis<Complex>>, 1> tmp = \
		array<shared_ptr<const ABasis<Complex>>, 1>({ std::move(bas) });
	discr0 = make_shared<AProjDiscr<Complex, 1, Complex>>(tmp);
	bas0 = discr0->bases[0].get();
	//throw away basis "bas" functions associated with last point
	nFree0 = bas0->getNCoef() - bas0->nFuncOnInterval()/2;
	assert(bas0->getBCTypeRight() == none);
	f0 = new Function<Complex, 1>(discr0);
	for (Int k = nFree0; k < bas0->getNCoef(); k++)
		f0->coef[k] = zzero;

//make additional functions part
	nBas1 = pn.size();
//	assert(nBas1 != 0);
	assert((cwf.size() == nBas1) && (y0nl.size() == nBas1));
	makeAddFunc(std::move(cwf), pn, y0nl);

//make basis
	nBas = nFree = nFree0 + nBas1;
	if (bas0->getBCTypeLeft() != none)
		nBas++; //NB! I assume that BC fixes 1 coefficient
	bcTypeLeft = bas0->getBCTypeLeft();
	bcTypeRight = none;
	bcLeft = bas0->getBCLeft();
	bcRight = zzero;
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->f(x, n);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] >= iInt)
			return rnl[n1](array<double, 1>({ x }));
		else //inl[n1] < iInt
			return fnl[n1].cwf->ulp(fnl[n1].pn*x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->d(x, n);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] >= iInt)
			return rnl[n1].df(array<double, 1>({ x }));
		else //inl[n1] < iInt
			return fnl[n1].pn*fnl[n1].cwf->dulp(fnl[n1].pn*x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->dd(x, n);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] >= iInt)
			return rnl[n1].d2f(array<double, 1>({ x }));
		else //inl[n1] < iInt
			return fnl[n1].pn*fnl[n1].pn*fnl[n1].cwf->ddulp(fnl[n1].pn*x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int k = 0; k < nFree0; k++)
		f0->coef[k] = func.coef[k];
	Complex res = bas0->f(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * fnl[n1].cwf->ulp(fnl[n1].pn*x);
		it++;
	}
	while (it != inl.end()) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * rnl[n1](array<double, 1>({ x }));
		it++;
	}

	return res;
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int k = 0; k < nFree0; k++)
		f0->coef[k] = func.coef[k];
	Complex res = bas0->d(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * fnl[n1].pn * fnl[n1].cwf->dulp(fnl[n1].pn*x);
		it++;
	}
	while (it != inl.end()) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * rnl[n1].df(array<double, 1>({ x }));
		it++;
	}

	return res;
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int k = 0; k < nFree0; k++)
		f0->coef[k] = func.coef[k];
	Complex res = bas0->dd(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * fnl[n1].pn * fnl[n1].pn * fnl[n1].cwf->ddulp(fnl[n1].pn*x);
		it++;
	}
	while (it != inl.end()) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * rnl[n1].d2f(array<double, 1>({ x }));
		it++;
	}

	return res;
}

template <class BasT>
void HybridBasis<BasT>::getNonzero(const double x, NumbersList &limits) const {
	if (x < grid.getLeftmostPoint() || x > grid.getRightmostPoint()) {
		limits.resize(0);
		return;
	}

	Int iInt = grid.getInterval(x); //number of interval enumerated from 0

	//"bas" part
	bas0->getNonzero(x, limits);
	if (iInt == (grid.getNIntervals() - 1)) { //last interval
		auto it = std::remove_if(limits.begin(), limits.end(), \
			[this](Int i) {return  i >= nFree0;});
		limits.erase(it, limits.end()); //C++20 replace with erase_if
	}

	//additional functions part
	for (Int k = nFree0; k < nFree; k++)
		limits.push_back(k);
}

template <class BasT>
Int HybridBasis<BasT>::nFuncOnInterval() const {
	assert(false);
	return 0;
}

template <>
Grid HybridBasis<AHermitSpline<Complex>>::collocGrid() const {
	const AHermitSpline<Complex> *bas0_ = \
		dynamic_cast<const AHermitSpline<Complex> *>(bas0);
	if (bas0_ == nullptr)
		assert(false);

	Int basPerNode = bas0_->nFuncOnInterval() / 2;
	Int iPoi = 0;
	Int nPoi = grid.getNIntervals()*basPerNode + nBas1;
	if (bas0_->getBCTypeLeft() != none) nPoi--;
	vector<double> poi; poi.resize(nPoi);

	Int nGaussP; double a, b, half, mid;
	vector<double> gaussP, gaussW;

	Int inl0;
	if (inl.size() > 0)
		inl0 = inl.front();
	else
		inl0 = 1;
	for (Int iInt = 0; iInt < grid.getNIntervals(); iInt++) {
		nGaussP = basPerNode;
		if (iInt == inl0 && bas0_->getBCTypeLeft() != none)
			nGaussP--; //bcTypeLeft = Dir or Neu or Mix
		nGaussP += countSorted(inl.begin(), inl.end(), iInt);
		gaussP.resize(nGaussP); gaussW.resize(nGaussP);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		a = grid[iInt]; b = grid[iInt + 1];
		half = 0.5*(b - a); mid = 0.5*(b + a);
		for (Int i = 0; i < gaussP.size(); i++)
			poi[iPoi++] = mid + half * gaussP[i];
		//Int cnt = countSorted(inl.begin(), inl.end(), iInt);
		//if (cnt > 0) {
		//	if (cnt == 1)
		//		poi[iPoi++] = b;
		//	else
		//		assert(false);
		//}
	}

	return Grid(poi);
}

template <class BasT>
GenMatrix<Complex> HybridBasis<BasT>::getOverlap() const {
	assert(false);
	GenMatrix<Complex> res(0);
	return res;
}

template <class BasT>
GenMatrix<Complex> HybridBasis<BasT>::getOverlapW(\
	std::function<double(double)> rho, const Int fdeg) const {
	assert(false);
	GenMatrix<Complex> res(0);
	return res;
}

template <class BasT>
void HybridBasis<BasT>::makeAddFunc(std::vector<unique_ptr<CWFCalculator>> &&cwf, \
	const std::vector<double> &pn, const std::vector<double> &y0nl) {
//find I_{nl}
	inl.reserve(nBas1);

	Int iInt;
	for (auto y0 : y0nl) {
		iInt = grid.getInterval(y0);
		assert(iInt != grid.getNIntervals()-1); //y0_{nl} cannot lie in the last interval
		inl.push_back(iInt);
	}

//sort functions with respect to I_{nl}
	std::vector<Int> p; p.resize(nBas1);
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(), \
			[this](Int i, Int j) {return inl[i] < inl[j]; });
	std::vector<Int> pinv; pinv.resize(nBas1);
	for (Int k = 0; k < nBas1; k++)
		pinv[p[k]] = k;

	fnl.reserve(nBas1);
	for (Int k = 0; k < nBas1; k++)
		fnl.push_back(Fnl(std::move(cwf[pinv[k]]), pn[pinv[k]]));
	std::sort(inl.begin(), inl.end());

//make R_{nl} functions
	makeRnl();
}

template <>
void HybridBasis<AHermitSpline<Complex>>::makeRnl() {
	const AHermitSpline<Complex> *bas0_ = \
		dynamic_cast<const AHermitSpline<Complex> *>(bas0);
	if (bas0_ == nullptr)
		assert(false);

	rnl.reserve(nBas1);
	Int nDer = bas0_->getNDeriv();
	assert(nDer <= 2); //only second derivative supported atm
	vector<Complex> val; val.resize(2*(nDer + 1));
	vector<double> poi; poi.resize(2);
	std::fill(val.begin(), val.begin()+nDer+1, zzero);
	poi[0] = grid.getLeftmostPoint();
	for (Int k = 0; k < nBas1; k++) {
		Grid gr(grid.getLeftmostPoint(), grid[inl[k]+1], 2);
		array<shared_ptr<const ABasis<Complex>>, 1> tmp;
		tmp[0] = make_shared<SimplePoly<Complex>>(gr, bas0_->getDeg());
		shared_ptr<AProjDiscr<Complex, 1, Complex>> discr = \
			make_shared<AProjDiscr<Complex, 1, Complex>>(tmp);
		rnl.push_back(Function<Complex, 1>(discr));
		poi[1] = grid[inl[k] + 1];
		for (Int iDer = 0; iDer <= nDer; iDer++) {
			switch (iDer) {
			case 0:
				val[nDer + 1] = fnl[k].cwf->ulp(fnl[k].pn*poi[1]);
				break;
			case 1:
				val[nDer + 2] = fnl[k].pn*fnl[k].cwf->dulp(fnl[k].pn*poi[1]);
				break;
			case 2:
				val[nDer + 3] = fnl[k].pn*fnl[k].pn*fnl[k].cwf->ddulp(fnl[k].pn*poi[1]);
				break;
			}
		}
		const SimplePoly<Complex> *sp = \
			dynamic_cast<const SimplePoly<Complex> *>(discr->bases[0].get());
		sp->interpH(rnl.back(), poi, val);
	}

}

template <class BasT>
HybridBasis<BasT>::~HybridBasis() {
	if (f0 != nullptr)
		delete f0;
}
*/

//var 2 --- Rnl support on 1 interval
/*
#include "HybridBasis.h"
#include "AHermitSpline.h"
#include "Algorithms.h"

template class HybridBasis<AHermitSpline<Complex>>;

template <class BasT>
HybridBasis<BasT>::HybridBasis(unique_ptr<const ABasis<Complex>> &&bas, \
	std::vector<unique_ptr<CWFCalculator>> &&cwf, \
	const std::vector<double> &pn, const std::vector<double> &y0nl) : \
		ABasis(bas->getGrid()) {
//grid
	assert(grid.getNIntervals() > 1); // >= 2 intervals are assumed
//make "bas" part
	array<unique_ptr<const ABasis<Complex>>, 1> tmp = \
		array<unique_ptr<const ABasis<Complex>>, 1>({ std::move(bas) });
	discr0 = make_shared<AProjDiscr<Complex, 1, Complex>>(tmp);
	bas0 = discr0->bases[0].get();
	//throw away basis "bas" functions associated with last point
	nFree0 = bas0->getNCoef() - bas0->nFuncOnInterval()/2;
	assert(bas0->getBCTypeRight() == none);
	f0 = new Function<Complex, 1>(discr0);
	for (Int k = nFree0; k < bas0->getNCoef(); k++)
		f0->coef[k] = zzero;

//make additional functions part
	nBas1 = pn.size();
	assert(nBas1 != 0);
	assert((cwf.size() == nBas1) && (y0nl.size() == nBas1));
	makeAddFunc(std::move(cwf), pn, y0nl);

//make basis
	nBas = nFree = nFree0 + nBas1;
	if (bas0->getBCTypeLeft() != none)
		nBas++; //NB! I assume that BC fixes 1 coefficient
	bcTypeLeft = bas0->getBCTypeLeft();
	bcTypeRight = none;
	bcLeft = bas0->getBCLeft();
	bcRight = zzero;
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->f(x, n);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] > iInt)
			return zzero;
		else if (inl[n1] == iInt)
			return rnl[n1](array<double, 1>({ x }))*fnl[n1].cwf->ulp(fnl[n1].pn*x);
		else //inl[n1] < iInt
			return fnl[n1].cwf->ulp(fnl[n1].pn*x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->d(x, n);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] > iInt)
			return zzero;
		else if (inl[n1] == iInt) {
			array<double, 1> xarr = array<double, 1>({ x });
			return rnl[n1](xarr)*fnl[n1].pn*fnl[n1].cwf->dulp(fnl[n1].pn*x) + \
				rnl[n1].df(xarr)*fnl[n1].cwf->ulp(fnl[n1].pn*x);
		}
		else //inl[n1] < iInt
			return fnl[n1].pn * fnl[n1].cwf->dulp(fnl[n1].pn*x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const Int n) const {
	if (n < nFree0) //"bas" function
		return bas0->dd(x, n);
	else { //additional function
		Int iInt = grid.getInterval(x); //number of interval enumerated from 0
		Int n1 = n - nFree0;
		if (inl[n1] > iInt)
			return zzero;
		else if (inl[n1] == iInt) {
			array<double, 1> xarr = array<double, 1>({ x });
			return rnl[n1](xarr)*fnl[n1].pn*fnl[n1].pn*fnl[n1].cwf->ddulp(fnl[n1].pn*x) + \
				2.0*rnl[n1].df(xarr)*fnl[n1].pn*fnl[n1].cwf->dulp(fnl[n1].pn*x) + \
				rnl[n1].d2f(xarr)*fnl[n1].cwf->ulp(fnl[n1].pn*x);
		}
		else //inl[n1] < iInt
			return fnl[n1].pn * fnl[n1].pn * fnl[n1].cwf->ddulp(fnl[n1].pn*x);
	}
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::f(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int k = 0; k < nFree0; k++)
		f0->coef[k] = func.coef[k];
	Complex res = bas0->f(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * fnl[n1].cwf->ulp(fnl[n1].pn*x);
		it++;
	}
	while (it != inl.end() && *it == iInt) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * \
			rnl[n1](array<double, 1>({ x }))*fnl[n1].cwf->ulp(fnl[n1].pn*x);
		it++;
	}

	return res;
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::d(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int k = 0; k < nFree0; k++)
		f0->coef[k] = func.coef[k];
	Complex res = bas0->d(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * fnl[n1].pn * fnl[n1].cwf->dulp(fnl[n1].pn*x);
		it++;
	}
	while (it != inl.end() && *it == iInt) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * \
			( rnl[n1](array<double, 1>({ x }))*fnl[n1].pn*fnl[n1].cwf->dulp(fnl[n1].pn*x) + \
				rnl[n1].df(array<double, 1>({ x }))*fnl[n1].cwf->ulp(fnl[n1].pn*x));
		it++;
	}

	return res;
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const AFunction<double> &func) const {
	assert(false);
}

template <class BasT>
Complex HybridBasis<BasT>::dd(const double x, const AFunction<Complex> &func) const {
	//"bas" part
	for (Int k = 0; k < nFree0; k++)
		f0->coef[k] = func.coef[k];
	Complex res = bas0->dd(x, *f0);

	//additional part
	Int iInt = grid.getInterval(x); //number of interval enumerated from 0
	auto lo = std::lower_bound(inl.begin(), inl.end(), iInt);
	Int n1;
	auto it = inl.begin();
	while (it != lo) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * fnl[n1].pn * fnl[n1].pn * fnl[n1].cwf->ddulp(fnl[n1].pn*x);
		it++;
	}
	while (it != inl.end() && *it == iInt) {
		n1 = it - inl.begin();
		res += func.coef[nFree0 + n1] * \
			(rnl[n1](array<double, 1>({ x }))*fnl[n1].pn*fnl[n1].pn*fnl[n1].cwf->ddulp(fnl[n1].pn*x) + \
				2.0*rnl[n1].df(array<double, 1>({ x }))*fnl[n1].pn*fnl[n1].cwf->dulp(fnl[n1].pn*x) + \
				rnl[n1].d2f(array<double, 1>({ x }))*fnl[n1].cwf->ulp(fnl[n1].pn*x));
		it++;
	}

	return res;
}

template <class BasT>
void HybridBasis<BasT>::getNonzero(const double x, NumbersList &limits) const {
	if (x < grid.getLeftmostPoint() || x > grid.getRightmostPoint()) {
		limits.resize(0);
		return;
	}

	Int iInt = grid.getInterval(x); //number of interval enumerated from 0

	//"bas" part
	bas0->getNonzero(x, limits);
	if (iInt == (grid.getNIntervals() - 1)) { //last interval
		auto it = std::remove_if(limits.begin(), limits.end(), \
			[this](Int i) {return  i >= nFree0;});
		limits.erase(it, limits.end()); //C++20 replace with erase_if
	}

	//additional functions part
	auto lb = std::upper_bound(inl.begin(), inl.end(), iInt);
	for (Int k = 0; k < lb-inl.begin(); k++)
		limits.push_back(nFree0+k);
}

template <class BasT>
Int HybridBasis<BasT>::nFuncOnInterval() const {
	assert(false);
	return 0;
}

template <>
Grid HybridBasis<AHermitSpline<Complex>>::collocGrid() const {
	const AHermitSpline<Complex> *bas0_ = \
		dynamic_cast<const AHermitSpline<Complex> *>(bas0);
	if (bas0_ == nullptr)
		assert(false);

	Int basPerNode = bas0_->nFuncOnInterval() / 2;
	Int iPoi = 0;
	Int nPoi = grid.getNIntervals()*basPerNode + nBas1;
	if (bas0_->getBCTypeLeft() != none) nPoi--;
	vector<double> poi; poi.resize(nPoi);

	Int nGaussP; double a, b, half, mid;
	vector<double> gaussP, gaussW;
	//first interval
	nGaussP = basPerNode;
	if (bas0_->getBCTypeLeft() != none) nGaussP--; //bcTypeLeft = Dir or Neu or Mix
	nGaussP += countSorted(inl.begin(), inl.end(), 0);
	gaussP.resize(nGaussP); gaussW.resize(nGaussP);
	gauleg(-1.0, 1.0, gaussP, gaussW);
	a = grid[0]; b = grid[1];
	half = 0.5*(b - a); mid = 0.5*(b + a);
	for (Int i = 0; i < gaussP.size(); i++)
		poi[iPoi++] = mid + half * gaussP[i];

	//inner and last intervals
	for (Int iInt = 1; iInt < grid.getNIntervals(); iInt++) {
		nGaussP = basPerNode;
		nGaussP += countSorted(inl.begin(), inl.end(), iInt);
		gaussP.resize(nGaussP); gaussW.resize(nGaussP);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		a = grid[iInt]; b = grid[iInt + 1];
		half = 0.5*(b - a); mid = 0.5*(b + a);
		for (Int i = 0; i < gaussP.size(); i++)
			poi[iPoi++] = mid + half * gaussP[i];
	}

	return Grid(poi);
}

template <class BasT>
GenMatrix<Complex> HybridBasis<BasT>::getOverlap() const {
	assert(false);
	GenMatrix<Complex> res(0);
	return res;
}

template <class BasT>
GenMatrix<Complex> HybridBasis<BasT>::getOverlapW(\
	std::function<double(double)> rho, const Int fdeg) const {
	assert(false);
	GenMatrix<Complex> res(0);
	return res;
}

template <class BasT>
void HybridBasis<BasT>::makeAddFunc(std::vector<unique_ptr<CWFCalculator>> &&cwf, \
	const std::vector<double> &pn, const std::vector<double> &y0nl) {
//find I_{nl}
	inl.reserve(nBas1);

	Int iInt;
	for (auto y0 : y0nl) {
		iInt = grid.getInterval(y0);
		assert(iInt != grid.getNIntervals()-1); //y0_{nl} cannot lie in the last interval
		inl.push_back(iInt);
	}

//sort functions with respect to I_{nl}
	std::vector<Int> p; p.resize(nBas1);
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(), \
			[this](Int i, Int j) {return inl[i] < inl[j]; });
	std::vector<Int> pinv; pinv.resize(nBas1);
	for (Int k = 0; k < nBas1; k++)
		pinv[p[k]] = k;

	fnl.reserve(nBas1);
	for (Int k = 0; k < nBas1; k++)
		fnl.push_back(Fnl(std::move(cwf[pinv[k]]), pn[pinv[k]]));
	std::sort(inl.begin(), inl.end());

//make R_{nl} functions
	makeRnl();
}

template <>
void HybridBasis<AHermitSpline<Complex>>::makeRnl() {
	const AHermitSpline<Complex> *bas0_ = \
		dynamic_cast<const AHermitSpline<Complex> *>(bas0);
	if (bas0_ == nullptr)
		assert(false);

	rnl.reserve(nBas1);
	for (Int k = 0; k < nBas1; k++)
		rnl.push_back(Function<Complex, 1>(discr0));

	Int nDer = bas0_->getNDeriv();
	Int iInt;
	vector<Complex> val; val.resize(grid.getNPoints()*(nDer+1));
	for (Int k = 0; k < rnl.size(); k++) {
		iInt = inl[k];
		std::fill(val.begin(), val.end(), zzero);
		for (Int i = iInt; i < grid.getNIntervals(); i++)
			val[(i+1)*(nDer+1)] = zone;

		bas0_->interpH(rnl[k], val);
	}

}

template <class BasT>
HybridBasis<BasT>::~HybridBasis() {
	if (f0 != nullptr)
		delete f0;
}
*/

