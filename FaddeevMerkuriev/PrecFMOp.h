#pragma once

#include "VectOperator.h"
#include "FMOp.h"

class PrecFMOp :
	public VectOperator<Complex, 3> {
public:
	unique_ptr<FMOp> fmop;
	unique_ptr<AnOperator> apprx;
	unique_ptr<AnOperator> prec;
	PrecFMOp(unique_ptr<FMOp> &&fmop, \
			unique_ptr<AnOperator> &&apprx, \
			unique_ptr<AnOperator> &&prec = nullptr) \
		: VectOperator<Complex, 3>(fmop->getDiscr()), \
			prec(std::move(prec)), fmop(std::move(fmop)), \
				apprx(std::move(apprx)) { }
	inline void times(const AFunction<Complex> &u, \
		AFunction<Complex> &res) override;
	inline void solve(AFunction<Complex> &rhssol) override;
	~PrecFMOp() = default;
};

inline void PrecFMOp::times(const AFunction<Complex> &u, \
	AFunction<Complex> &res) {
	fmop->times(u, res);
	if (prec == nullptr)
		apprx->solve(res);
	else {
		unique_ptr<AFunction<Complex>> tmp \
			= res.clone();
		prec->times(*tmp, res);
	}
}

inline void PrecFMOp::solve(AFunction<Complex> &rhssol) {
	assert(false);
}