#include "Particle.h"

Particle::Particle(const double &mass, const double &charge) \
			: m(mass), q(charge) { }

void Particle::print() const {
	cout << "Particle:  mass = " << m << ",  charge = " << q;
	cout << endl;
}

bool Particle::operator==(const Particle &rhs) const {
	return (m == rhs.m) && (q == rhs.q);
}

Particle::~Particle(void) { }
