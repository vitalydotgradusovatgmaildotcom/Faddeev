#pragma once

#include "Particle.h"
#include "CoulombPot.h"

//2 body subsystem
class Pair {
public:
	Int alpha; //pair number
	std::vector<Particle> p; //Particles
	//short-range potential
	shared_ptr<APotential<Complex>> vs = nullptr;
	//coulomb potential
	shared_ptr<CoulombPot<Complex>> vc = nullptr;
	double mu; //reduced mass
	double q12;
	bool id; //true if particles are identical
	Pair(const Pair &rhs) = delete;
	Pair & operator=(const Pair &rhs) = delete;
	Pair & operator=(Pair &&rhs) = delete;
	//physical 2 reduced coordinate
	void reduceCoo(double &x) const {
		x *= sqrt(2.0*mu);
	};
	~Pair() = default;
protected:
	Pair() = default;
	Pair(Pair &&rhs) = default;
	friend class System3Body;
};

