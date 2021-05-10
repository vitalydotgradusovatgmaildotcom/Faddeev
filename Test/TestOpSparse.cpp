#include "TestOpSparse.h"


TestOpSparse::TestOpSparse(const shared_ptr<const ADiscretization<TEST_OP_SP_TYPE, 1>> &discr) \
				: Operator<TEST_OP_SP_TYPE, 1>(discr), \
					matr(discr->getN(), 100*discr->getN()) { }

Function<TEST_OP_SP_TYPE, 1> TestOpSparse::operator*(const Function<TEST_OP_SP_TYPE, 1> &u) const {
	Function<TEST_OP_SP_TYPE, 1> res(u);
	res.coef = matr*u.coef;
	return res;
}

void TestOpSparse::times(const AFunction<TEST_OP_SP_TYPE> &u, AFunction<TEST_OP_SP_TYPE> &res) {
	const Function<TEST_OP_SP_TYPE, 1> & u_ = \
		dynamic_cast<const Function<TEST_OP_SP_TYPE, 1> &>(u);
	Function<TEST_OP_SP_TYPE, 1> & res_ = \
		dynamic_cast<Function<TEST_OP_SP_TYPE, 1> &>(res);
	res_ = *this * u_;
}

void TestOpSparse::solve(AFunction<Complex> &rhssol) {
	assert(false);
}

TestOpSparse::~TestOpSparse(void) { }
