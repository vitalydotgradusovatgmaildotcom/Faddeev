#pragma once

#include "VectFunction.h"
#include "CompsDiscr.h"

class Components :
	public VectFunction<Complex, 3> {
public:
	Components(const shared_ptr<const CompsDiscr> &cd);
	Components(const Components &rhs) = default;
	Components(Components &&rhs) = default;
	Components & operator=(const Components &rhs) = default;
	Components & operator=(const AFunction<Complex> &rhs) override;
	Components & operator=(Components &&rhs) = default;
	Components & operator=(AFunction<Complex> &&rhs) override;
	Components operator*(const Complex x) const;
	Components operator+(const Components &u) const;
	Components operator-(const Components &u) const;
	inline const Function<Complex, 3> & getF(const Int alpha, const Int M) const;
	unique_ptr<AFunction<Complex>> clone() const override;
	//shared_ptr<const CompsDiscr> getDiscr() const;
	~Components(void);
protected:
	shared_ptr<const CompsDiscr> cd;
};

inline const Function<Complex, 3> & Components::getF(const Int alpha, const Int M) const {
	return VectFunction<Complex, 3>::getF(cd->alphaM2IFunc(alpha, M));
}

