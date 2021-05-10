#include "TestOperator.h"


TestOperator::TestOperator(const shared_ptr<const ADiscretization<TEST_OP_TYPE, 1>> &discr) \
				: Operator<TEST_OP_TYPE, 1>(discr) {
	matr.resize(discr->getN());
}

/*
TestOperator & TestOperator::operator=(const AnOperator<TEST_OP_TYPE> &rhs) {
	Operator::operator=(rhs);
	const TestOperator & rhs_ = dynamic_cast<const TestOperator &>(rhs);
	this->matr = rhs_.matr;

	return *this;
}
*/

Function<TEST_OP_TYPE, 1> TestOperator::operator*(const Function<TEST_OP_TYPE, 1> &u) const {
	Function<TEST_OP_TYPE, 1> res(u);
	res.coef = matr*u.coef;
	return res;
}

void TestOperator::times(const AFunction<TEST_OP_TYPE> &u, AFunction<TEST_OP_TYPE> &res) {
	const Function<TEST_OP_TYPE, 1> & u_ = \
		dynamic_cast<const Function<TEST_OP_TYPE, 1> &>(u);
	Function<TEST_OP_TYPE, 1> & res_ = \
		dynamic_cast<Function<TEST_OP_TYPE, 1> &>(res);
	res_ = *this * u_;
}

void TestOperator::solve(AFunction<Complex> &rhssol) {
	assert(false);
}

TestOperator::~TestOperator(void) { }
