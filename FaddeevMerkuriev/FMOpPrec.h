#pragma once

#include "VectOperator.h"
#include "FMOp.h"

class FMOpPrec :
	public VectOperator<Complex, 3> {
public:
	FMOpPrec(const FMOp &fmop) : \
		VectOperator<Complex, 3>(fmop.getDiscr()), \
		matr(ILU0(*fmop.matr)) { };
	inline void times(const AFunction<Complex> &u, \
		AFunction<Complex> &res) override;
	inline void solve(AFunction<Complex> &rhssol) override;
	~FMOpPrec() = default;
protected:
	SparseMatr<Complex> matr;
};

inline void FMOpPrec::times(const AFunction<Complex> &u, \
	AFunction<Complex> &res) {
	res = u;
	matr.solve(res.coef);
}

inline void FMOpPrec::solve(AFunction<Complex> &rhssol) {
	assert(false);
}