#pragma once

#include "VectDiscr.h"

class Components;

class CompsDiscr :
	public VectDiscr<Complex, 3> {
public:
	CompsDiscr(const vector<shared_ptr<const ADiscretization<Complex, 3>>> &discrs, \
				const Int nComp, const Int Mmin, const Int Mmax, \
					const vector<Int> &solComps);
	inline Int getRaw(const Int alpha, const Int M, const std::array<Int, 3> &is) const;
	inline void getInds(const Int i, Int &alpha, Int &M, std::array<Int, 3> &is) const;
	inline Int getNComp() const;
	vector<Int> getSolComps() const; //solution components numbers
	inline const ADiscretization<Complex, 3> & get(const Int alpha, const Int M) const;
	~CompsDiscr(void) = default;
	friend class Components;
protected:
	Int nComp;
	vector<Int> solComps; //solution components
	vector<Int> comp2sol; //components to solution components
	Int Mmin, Mmax, nM;
	void makeMap();
	inline Int alphaM2IFunc(const Int alpha, const Int M) const;
	inline void iFunc2AlphaM(const Int iFunc, Int &alpha, Int &M) const;
};

inline Int CompsDiscr::getRaw(const Int alpha, const Int M, const std::array<Int, 3> &is) const {
	return VectDiscr<Complex, 3>::getRaw(alphaM2IFunc(comp2sol[alpha], M), is);
}

//NB! returns SOME component alpha that maps to solution component
inline void CompsDiscr::getInds(const Int i, Int &alpha, Int &M, std::array<Int, 3> &is) const {
	static Int iFunc;
	VectDiscr<Complex, 3>::getInds(i, iFunc, is);
	iFunc2AlphaM(iFunc, alpha, M);
	alpha = solComps[alpha];
}

inline Int CompsDiscr::alphaM2IFunc(const Int alpha, const Int M) const {
	return alpha*nM + M - Mmin;
}

inline void CompsDiscr::iFunc2AlphaM(const Int iFunc, Int &alpha, Int &M) const {
	static auto qr = std::div(iFunc, nM);
	alpha = qr.quot;
	M = qr.rem + Mmin;
}

inline Int CompsDiscr::getNComp() const {
	return nComp;
}

inline const ADiscretization<Complex, 3> & CompsDiscr::get(const Int alpha, const Int M) const {
	return (*this)[alphaM2IFunc(comp2sol[alpha], M)];
}