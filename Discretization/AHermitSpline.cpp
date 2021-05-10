#include "AHermitSpline.h"

template class AHermitSpline<double>;
template class AHermitSpline<Complex>;

template <class TBas>
AHermitSpline<TBas>::AHermitSpline(const Grid &gr) \
		: basPerNode(0), ABasis<TBas>(gr) { }

template <class TBas>
AHermitSpline<TBas>::AHermitSpline(const Grid &gr, \
		const enum bc bcTypeLeft, const enum bc bcTypeRight, \
		const TBas bcLeft, const TBas bcRight) \
		: basPerNode(0), ABasis<TBas>(gr, bcTypeLeft, bcTypeRight, bcLeft, bcRight) {}

template <class TBas>
void AHermitSpline<TBas>::getNonzero(const double x, NumbersList &limits) const {

	if (x < grid.getLeftmostPoint() || x > grid.getRightmostPoint()) {
		limits.resize(0);
		return;
	}

	Int iInt = grid.getInterval(x); //number of Interval enumerated from 0
	Int iLeft = basPerNode * iInt; //enumerate from 0

	limits.resize(0);
	Int sh = 0;
	if (bcTypeLeft != none)
		sh++;
	if (iInt == 0) {
		if (grid.getNIntervals() == 1 && bcTypeRight != none)
			sh++;
		for (Int k = 0; k < nFuncOnInterval()-sh; k++)
			limits.push_back(k);
	} else if ( iInt == (grid.getNIntervals()-1) ) {
		Int k = iLeft - sh;
		if (bcTypeRight != none)
			sh++;
		for (; k < nBas-sh; k++)
			limits.push_back(k);
	} else { //Internal Interval
		for (Int i = 0; i < nFuncOnInterval(); i++)
			limits.push_back(iLeft + i - sh);
	}
}

template <class TBas>
void AHermitSpline<TBas>::initRefPoints() {

	lPoi.reserve(nBas); mPoi.reserve(nBas); rPoi.reserve(nBas);

	Int iNum; //Interval number 0, 1, ...
	for (Int k = 0; k < basPerNode; k++) {
		lPoi.push_back(grid[0]);
		mPoi.push_back(grid[0]);
		rPoi.push_back(grid[1]);
	}
	for (Int k = basPerNode; k < nBas-basPerNode; k++) {
		iNum = k / basPerNode;
		lPoi.push_back(grid[iNum-1]);
		mPoi.push_back(grid[iNum]);
		rPoi.push_back(grid[iNum+1]);
	}
	for (Int k = nBas-basPerNode; k < nBas; k++) {
		lPoi.push_back(grid[grid.getNPoints()-2]);
		mPoi.push_back(grid[grid.getNPoints()-1]);
		rPoi.push_back(grid[grid.getNPoints()-1]);
	}
}

template <class TBas>
Grid AHermitSpline<TBas>::collocGrid() const {
	Int iPoi = 0;
	Int nPoi = basPerNode*(grid.getNIntervals()+1);
	if (bcTypeLeft != none) nPoi--;
	if (bcTypeRight != none) nPoi--;
	vector<double> poi; poi.resize(nPoi);

	Int nGaussP; double a, b, half, mid;
	vector<double> gaussP, gaussW;
//1 interval
	if (grid.getNIntervals() == 1) {
		nGaussP = 2*basPerNode;
		if (bcTypeLeft != none) nGaussP--;
		if (bcTypeRight != none) nGaussP--;
		gaussP.resize(nGaussP); gaussW.resize(nGaussP);
		gauleg(-1.0, 1.0, gaussP, gaussW);
		a = grid[0]; b = grid[1];
		half = 0.5*(b-a); mid = 0.5*(b+a);
		for (Int i = 0; i < gaussP.size(); i++)
			poi[iPoi++] = mid + half*gaussP[i];
		return Grid(poi);
	} //end 1 interval

// > 1 intervals
	//first interval
	nGaussP = basPerNode+(basPerNode+1)/2;
	if (bcTypeLeft != none) nGaussP--; //bcTypeLeft = Dir or Neu or Mix
	gaussP.resize(nGaussP); gaussW.resize(nGaussP);
	gauleg(-1.0, 1.0, gaussP, gaussW);
	a = grid[0]; b = grid[1];
	half = 0.5*(b-a); mid = 0.5*(b+a);
	for (Int i = 0; i < gaussP.size(); i++)
		poi[iPoi++] = mid + half*gaussP[i];

	//inner intervals
	nGaussP = basPerNode;
	gaussP.resize(nGaussP); gaussW.resize(nGaussP);
	gauleg(-1.0, 1.0, gaussP, gaussW);
	for (Int iInt = 1; iInt < grid.getNIntervals()-1; iInt++) {
		a = grid[iInt]; b = grid[iInt+1];
		half = 0.5*(b-a); mid = 0.5*(b+a);
		for (Int i = 0; i < gaussP.size(); i++)
			poi[iPoi++] = mid + half*gaussP[i];
	}

	//last interval
	nGaussP = basPerNode+basPerNode/2;
	if (bcTypeRight != none) nGaussP--; //bcTypeRight = Dir or Neu or Mix
	gaussP.resize(nGaussP); gaussW.resize(nGaussP);
	gauleg(-1.0, 1.0, gaussP, gaussW);
	a = grid[grid.getNPoints()-2]; b = grid[grid.getNPoints()-1];
	half = 0.5*(b-a); mid = 0.5*(b+a);
	for (Int i = 0; i < gaussP.size(); i++)
		poi[iPoi++] = mid + half*gaussP[i];

	return Grid(poi);
}

template <class TBas>
GenMatrix<TBas> AHermitSpline<TBas>::getOverlap() const {
	GenMatrix<TBas> ovlpM(nFree);
	ovlpM.fill(TBas());

	vector<double> gaussX(deg+1), gaussW(deg+1);
	gauleg(-1.0, 1.0, gaussX, gaussW);

	double a, b, half, mid, t;
	NumbersList limits_;
	for (Int iInt = 0; iInt < grid.getNIntervals(); iInt++) {
		a = grid[iInt]; b = grid[iInt+1];
		half = 0.5*(b-a); mid = 0.5*(b+a);
		getNonzero(mid, limits_);
		for (Int k = 0; k < limits_.size(); k++)
			for (Int l = k; l < limits_.size(); l++) {
				for (Int i = 0; i < gaussX.size(); i++) {
					t = mid + half*gaussX[i];
					ovlpM[limits_[k]][limits_[l]] += \
						half * gaussW[i] * CONJ( this->f(t, limits_[k]) )*this->f(t, limits_[l]);
				}
				ovlpM[limits_[l]][limits_[k]] = \
						CONJ( ovlpM[limits_[k]][limits_[l]] );
			}
	}

	return ovlpM;
}

template <class TBas>
GenMatrix<TBas> AHermitSpline<TBas>::getOverlapW(\
	std::function<double(double)> rho, const Int fdeg) const {
	GenMatrix<TBas> ovlpM(nFree);
	ovlpM.fill(TBas());

	vector<double> gaussX(deg + 1 + (fdeg+1)/2), \
		gaussW(deg + 1 + (fdeg + 1) / 2);
	gauleg(-1.0, 1.0, gaussX, gaussW);

	double a, b, half, mid, t;
	NumbersList limits_;
	for (Int iInt = 0; iInt < grid.getNIntervals(); iInt++) {
		a = grid[iInt]; b = grid[iInt + 1];
		half = 0.5*(b - a); mid = 0.5*(b + a);
		getNonzero(mid, limits_);
		for (Int k = 0; k < limits_.size(); k++)
			for (Int l = k; l < limits_.size(); l++) {
				for (Int i = 0; i < gaussX.size(); i++) {
					t = mid + half * gaussX[i];
					ovlpM[limits_[k]][limits_[l]] += \
						half * gaussW[i] * CONJ( this->f(t, limits_[k]) )*this->f(t, limits_[l])*rho(t);
				}
				ovlpM[limits_[l]][limits_[k]] = \
					CONJ( ovlpM[limits_[k]][limits_[l]] );
			}
	}

	return ovlpM;
}

template <class TBas>
void AHermitSpline<TBas>::interpH(AFunction<TBas> &func, \
	const vector<TBas> &val) const {
	assert(val.size() == nBas);
	assert(func.size() == nFree);

	switch (bcTypeLeft) {
	case Dir:
		assert(val[0] == bcLeft);
		break;
	case Neu:
		assert(val[1] == bcLeft);
		break;
	case Mix:
		assert(val[1] / val[0] == bcLeft);
		break;
	}

	switch (bcTypeRight) {
	case Dir:
		assert(val[nBas-basPerNode] == bcRight);
		break;
	case Neu:
		assert(val[nBas-basPerNode+1] == bcRight);
		break;
	case Mix:
		assert(val[nBas-basPerNode+1] / val[nBas-basPerNode] == bcRight);
		break;
	}

	for (Int k = 0; k < nFree; k++)
		func[k] = val[toInternalN(k)];
}