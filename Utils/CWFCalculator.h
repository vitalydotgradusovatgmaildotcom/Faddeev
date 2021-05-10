#pragma once

#include "cwfcomp.H"

//Function definitions are those of Messiah
inline Complex sigmaL(const Complex l, const Complex eta) {
	return sigma_l_calc(l, eta);
}

class CWFCalculator {
public:
	CWFCalculator(const Complex l, const Complex eta) : \
		l(l), eta(eta), cwf(true, l, eta), eisigmal(exp(zi*sigmaL(l, eta))), \
			llp1(l*(l+1)), eta2(2.0*eta) { };
	CWFCalculator(const CWFCalculator &rhs) : l(rhs.l), eta(rhs.eta), \
		cwf(true, rhs.l, rhs.eta), eisigmal(rhs.eisigmal), \
			llp1(rhs.llp1), eta2(rhs.eta2) { };
	CWFCalculator(CWFCalculator &&rhs) : l(rhs.l), eta(rhs.eta), \
		cwf(true, rhs.l, rhs.eta), eisigmal(rhs.eisigmal), \
			llp1(rhs.llp1), eta2(rhs.eta2) { };
	CWFCalculator & operator=(const CWFCalculator &rhs) = delete;
	CWFCalculator & operator=(CWFCalculator &&rhs) = delete;
	Complex fl(Complex z) {
		Complex f, df;
		cwf.F_dF(z, f, df);
		return f;
	};
	Complex flExp(Complex z) {
		Complex f, df;
		cwf.F_dF(z, f, df);
		return eisigmal * f;
	};
	Complex gl(Complex z) {
		Complex f, df;
		cwf.G_dG(z, f, df);
		return f;
	};
	Complex ulp(Complex z) {
		Complex f, df;
		cwf.H_dH(1, z, f, df);
		return f / eisigmal;
	};
	Complex ulm(Complex z) {
		Complex f, df;
		cwf.H_dH(-1, z, f, df);
		return eisigmal * f;
	};
	Complex dfl(Complex z) {
		Complex f, df;
		cwf.F_dF(z, f, df);
		return df;
	};
	Complex dgl(Complex z) {
		Complex f, df;
		cwf.G_dG(z, f, df);
		return df;
	};
	Complex dulp(Complex z) {
		Complex f, df;
		cwf.H_dH(1, z, f, df);
		return df / eisigmal;
	};
	Complex dulm(Complex z) {
		Complex f, df;
		cwf.H_dH(-1, z, f, df);
		return eisigmal * df;
	};
	Complex ddfl(Complex z) {
		Complex f, df;
		cwf.F_dF(z, f, df);
		if (z == zzero) {
			assert(l.imag() == 0.0);
			Complex c0 = eta == zzero ? zone : sqrt(2*PI*eta/(exp(2*PI*eta)-1));
			if (l.real() == 0.0)
				return eta2 * c0;
			else if (l.real() == 1.0)
				return 2.0*c0*sqrt(1.0 + eta2) / 3;
			else
				return zzero;
		}
		return llp1*((f/z)/z) + eta2*(f/z) - f;
	};
	Complex ddgl(Complex z) {
		Complex f, df;
		cwf.G_dG(z, f, df);
		return llp1 * ((f / z) / z) + eta2 * (f / z) - f;
	};
	Complex ddulp(Complex z) {
		Complex f, df;
		cwf.H_dH(1, z, f, df);
		return (llp1 * ((f / z) / z) + eta2 * (f / z) - f) / eisigmal;
	};
	Complex ddulm(Complex z) {
		Complex f, df;
		cwf.H_dH(-1, z, f, df);
		return eisigmal * (llp1 * ((f / z) / z) + eta2 * (f / z) - f);
	};
	~CWFCalculator() = default;
protected:
	Coulomb_wave_functions cwf;
	Complex l, eta;
	Complex eisigmal;
	Complex llp1, eta2;
};

