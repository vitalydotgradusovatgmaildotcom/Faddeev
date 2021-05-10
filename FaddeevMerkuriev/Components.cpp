#include "Components.h"


Components::Components(const shared_ptr<const CompsDiscr> &cd) \
		: cd(cd), VectFunction<Complex, 3>(cd) { }

Components & Components::operator=(const AFunction<Complex> &rhs) {
	const Components &rhs_ = \
		dynamic_cast<const Components &>(rhs);
	this->operator=(rhs_);
	return *this;
}
/*
Components & Components::operator=(Components &&rhs) {
	VectFunction<Complex, 3>::operator=(std::move(rhs));
	cd = std::move(rhs.cd); rhs.cd.reset();
	return *this;
}*/

Components & Components::operator=(AFunction<Complex> &&rhs) {
	Components &&rhs_ = \
		dynamic_cast<Components &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

Components Components::operator*(const Complex x) const {
	Components res(*this);
	res *= x;
	return res;
}

Components Components::operator+(const Components &u) const {
	Components res(*this);
	res += u;
	return res;
}

Components Components::operator-(const Components &u) const {
	Components res(*this);
	res -= u;
	return res;
}

unique_ptr<AFunction<Complex>> Components::clone() const {
	return make_unique<Components>(*this);
}

/*
shared_ptr<const CompsDiscr> Components::getDiscr() const {
	return cd;
}*/

Components::~Components(void) { }
