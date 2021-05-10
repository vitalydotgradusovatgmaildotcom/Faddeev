#include "System3Body.h"
#include "HeTTY.h"
#include "KTTY.h"
#include "CoulSharplyCut.h"
#include "ExpPot.h"

System3Body::System3Body( \
	const string &name, const array<double, 3> &m, \
	const array<double, 3> &q, const string &idstr, \
		const Int J, const Int Mproj, const Int tau) \
		: system(name), J(J), Mproj(Mproj), tau(tau) {
	makeParticles(m, q);
	makeIdentical(idstr);
	calculateReducedMasses();
	calculateRotationMatrix();
	calculateCoulombCouplings();
	makePotentials();
}

void System3Body::makeParticles(const array<double, 3> &m, \
							const array<double, 3> &q) {
	p.reserve(3);
	for (Int k = 0; k < 3; k++)
		p.push_back(Particle(m[k], q[k]));
}

void System3Body::makeIdentical(const string &idstr) {
	if (idstr == "3 symmetric") {
		id = three_sym;
		assert(p[0] == p[1] && p[1] == p[2]);
	} else if (idstr == "3 antisymmetric") {
		id = three_asym;
		assert(p[0] == p[1] && p[1] == p[2]);
	} else if (idstr == "2 symmetric") {
		id = two_sym;
		assert(p[0] == p[1]);
	} else if (idstr == "2 antisymmetric") {
		id = two_asym;
		assert(p[0] == p[1]);
	} else
		id = no_ident;
}

void System3Body::calculateReducedMasses() {
	double m1, m2, m3, M;
	M = p[0].m + p[1].m + p[2].m;

    m1 = p[1].m; m2 = p[2].m; m3 = p[0].m;
    m12[0] = m1*m2/(m1+m2);
	m123[0] = m3*(m1 + m2)/M;
    //
    m1 = p[0].m; m2 = p[2].m; m3 = p[1].m;
    m12[1] = m1*m2/(m1+m2);
	m123[1] = m3*(m1 + m2)/M;
    //
    m1 = p[0].m; m2 = p[1].m; m3 = p[2].m;
    m12[2] = m1*m2/(m1+m2);
	m123[2] = m3*(m1 + m2)/M;
}

void System3Body::calculateRotationMatrix() {
	double M = p[0].m + p[1].m + p[2].m;
	for (Int alpha = 0; alpha < 3; alpha++)
		for (Int beta = 0; beta < 3; beta++) {
			c[beta][alpha] = -sqrt( p[beta].m*p[alpha].m/( (M-p[beta].m)*(M-p[alpha].m)) );
			//s[beta][alpha] = pow(-1.0, beta-alpha)*SGN(alpha-beta)*sqrt(1-c[beta][alpha]*c[beta][alpha]);
			s[beta][alpha] = pow(-1.0, beta - alpha)*SGN(beta - alpha)*sqrt(1 - c[beta][alpha] * c[beta][alpha]);
		}
	c[0][0] = 1.0; c[1][1] = 1.0; c[2][2] = 1.0;
	s[0][0] = 0.0; s[1][1] = 0.0; s[2][2] = 0.0;
	for (Int alpha = 0; alpha < 3; alpha++)
		for (Int beta = 0; beta < 3; beta++)
			c2_s2[beta][alpha] = c[beta][alpha]*c[beta][alpha]-s[beta][alpha]*s[beta][alpha];
}

void System3Body::calculateCoulombCouplings() {
	q12[0] = p[1].q*p[2].q;
	q12[1] = p[2].q*p[0].q;
	q12[2] = p[0].q*p[1].q;
	q123[0] = p[0].q*(p[1].q + p[2].q);
	q123[1] = p[1].q*(p[2].q + p[0].q);
	q123[2] = p[2].q*(p[0].q + p[1].q);
}

void System3Body::makePotentials() {
	//NB scaled Coulomb potentials!
	for (Int alpha = 0; alpha < 3; alpha++) {
		vc[alpha] = make_shared<CoulombPot<Complex>>(q12[alpha]);
		vs[alpha] = make_shared<ZeroPot<Complex>>();
	}
	//short-range
	if (system == "HeMol_TTY")
		for (Int alpha = 0; alpha < 3; alpha++)
			vs[alpha] = make_shared<HeTTY<Complex>>();
	else if (system == "He2Li_KTTY") {
		for (Int alpha = 0; alpha < 2; alpha++)
			vs[alpha] = make_shared<KTTY<Complex>>( \
				2.430857, 1.04911, 0.00381298, 22.507, 1083.2, 72602.1);
		vs[2] = make_shared<HeTTY<Complex>>();
	} else if (system == "He2Na_KTTY") {
		for (Int alpha = 0; alpha < 2; alpha++)
			vs[alpha] = make_shared<KTTY<Complex>>(\
				2.218564, 1.00872, 0.00399053, 23.768, 1307.6, 94563.2);
		vs[2] = make_shared<HeTTY<Complex>>();
	} else if (system == "He2K_KTTY") {
		for (Int alpha = 0; alpha < 2; alpha++)
			vs[alpha] = make_shared<KTTY<Complex>>(\
				1.568281, 0.86941, 0.00466213, 34.038, 2525.2, 237538.0);
		vs[2] = make_shared<HeTTY<Complex>>();
	} else if (system == "He2Rb_KTTY") {
		for (Int alpha = 0; alpha < 2; alpha++)
			vs[alpha] = make_shared<KTTY<Complex>>(\
				1.440646, 0.83839, 0.00482456, 36.289, 2979.0, 300406.0);
		vs[2] = make_shared<HeTTY<Complex>>();
	} else if (system == "He2Cs_KTTY") {
		for (Int alpha = 0; alpha < 2; alpha++)
			vs[alpha] = make_shared<KTTY<Complex>>(\
				1.440646, 0.83839, 0.00482456, 41.417, 3903.4, 453443.0);
		vs[2] = make_shared<HeTTY<Complex>>();
	} else if (system == "TestScatter") {
		for (Int alpha = 0; alpha < 3; alpha++)
			//vs[alpha] = make_shared<CoulSharplyCut<Complex>>(-q12[alpha], 2.5);
			vs[alpha] = make_shared<ExpPot<Complex>>(-100.0, 1.0/2.5);
	}

}

void System3Body::reduceCoo(double &x, double &y, const Int alpha) const {
	x *= sqrt(2.0*m12[alpha]);
	y *= sqrt(2.0*m123[alpha]);
}

void System3Body::rotReduced(const jacobiCoo &coo, const Int alpha, jacobiCoo &res, const Int beta) const {
	double o11x = c[beta][alpha]*coo.x;
	double o21x = -s[beta][alpha]*coo.x;
	double o12y = s[beta][alpha]*coo.y;
	double o22y = c[beta][alpha]*coo.y;

	res.x = sqrt(o11x*(o11x+2.0*o12y*coo.z)+o12y*o12y);
    res.y = sqrt(o21x*(o21x+2.0*o22y*coo.z)+o22y*o22y);
    res.z = (o11x*o21x+c2_s2[beta][alpha]*coo.x*coo.y*coo.z+ \
                        o12y*o22y)/(res.x*res.y);
}

Pair System3Body::getPair(const Int alpha) const {
	Pair res;

	
	res.id = false;
	if (id == three_sym || id == three_asym) {
		res.id = true;
	}
	if ((id == two_sym || id == two_asym) && alpha == 2) {
		res.id = true;
	}

	res.alpha = alpha;
	res.p.reserve(2);
	for (Int k = 0; k < 3; k++)
		if (k != alpha)
			res.p.push_back(p[k]);
	if (res.id)
		assert(res.p[0] == res.p[1]);
	res.mu = m12[alpha];
	res.q12 = q12[alpha];
	res.vc = vc[alpha];
	res.vs = vs[alpha];

	return res;
}

void System3Body::printSystem() const {
	cout << "System: " << system << endl << endl;
	
	for (Int k = 0; k < p.size(); k++)
		p[k].print();
	cout << endl;

	cout << "Reduced masses:" << endl;
	cout << "m_{12}: " << m12[0] << "  " << \
		m12[1] << "  " << m12[2] << "  " << endl;
	cout << "m_{123}: " << m123[0] << "  " << \
		m123[1] << "  " << m123[2] << "  " << endl;
	cout << endl;

	switch (id) {
	case three_sym:
		cout << "1, 2, 3 identical, symmetric w.f." << endl;
		break;
	case three_asym:
		cout << "1, 2, 3 identical, antisymmetric w.f." << endl;
		break;
	case two_sym:
		cout << "1 and 2 identical, symmetric w.f." << endl;
		break;
	case two_asym:
		cout << "1 and 2 identical, antisymmetric w.f." << endl;
		break;
	case no_ident:
		cout << "No identical" << endl;
		break;
	default:
		assert(false);
		break;
	}
	cout << endl;

	cout << "J = " << J << ", M = " << Mproj << \
		", tau (physical) = " << tau*pow(-1.0, J) << endl << endl;
}