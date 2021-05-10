#pragma once

#include "Particle.h"
#include "CoulombPot.h"
#include "Pair.h"

struct jacobiCoo {
	double x;
	double y;
	double z;
};

enum identical {no_ident, two_sym, two_asym, three_sym, three_asym};

class System3Body {
public:
	string system;
	std::vector<Particle> p; //particles
	//short-range potential
	std::array<shared_ptr<APotential<Complex>>, 3> vs \
		= {nullptr, nullptr, nullptr};
	//coulomb potential
	std::array<shared_ptr<CoulombPot<Complex>>, 3> vc \
		= { nullptr, nullptr, nullptr };
	std::array<double, 3> m12, m123;
	std::array<double, 3> q12, q123;
	enum identical id;
	Int J, Mproj, tau;
	double c2_s2[3][3], s[3][3], c[3][3];
	System3Body(const string &name, const array<double, 3> &m, \
		const array<double, 3> &q, const string &idstr, \
		const Int J, const Int Mproj, const Int tau);
	System3Body(const System3Body &rhs) = delete;
	System3Body(System3Body &&rhs) = delete;
	System3Body & operator=(const System3Body &rhs) = delete;
	System3Body & operator=(System3Body &&rhs) = delete;
	//physical 2 reduced Jacobi coordinates
	void reduceCoo(double &x, double &y, const Int alpha) const;
	//rotates reduced Jacobi coordinates
	void rotReduced(const jacobiCoo &coo, const Int alpha, jacobiCoo &res, const Int beta) const;
	Pair getPair(const Int alpha) const;
	void printSystem() const;
	~System3Body(void) = default;

protected:
	void makeParticles(const array<double, 3> &m, const array<double, 3> &q);
	void makeIdentical(const string &idstr);
	void calculateReducedMasses();
	void calculateRotationMatrix();
	void calculateCoulombCouplings();
	void makePotentials();
};

