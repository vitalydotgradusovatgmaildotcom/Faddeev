#pragma once

#include "cph.h"

class Particle {
public:
	double m;
	double q;
	Particle(const double &mass, const double &charge);
	Particle(const Particle &p) = default;
	Particle(Particle &&rhs) = default;
	Particle & operator=(const Particle &rhs) = default;
	Particle & operator=(Particle &&rhs) = default;
	bool operator==(const Particle &rhs) const;
	void print() const;
	~Particle(void);
};

