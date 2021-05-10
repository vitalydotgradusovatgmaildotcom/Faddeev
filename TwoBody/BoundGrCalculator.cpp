#include "BoundGrCalculator.h"
#include "MathAn.h"
#include "Mappings.h"


BoundGrCalculator::BoundGrCalculator(const Pair &pair) : \
		config(FaddeevConfigurator::getInstance()), pair(pair) {
	//energy conversion factor
	eConv = config.eConv;
	makeDiscretization();
	makePotentials();
	makeChi();
}

void BoundGrCalculator::makeDiscretization() {
	double xmax = config.xmax[pair.alpha];
	scaleCoo(xmax);
	gr = Grid(0.0, xmax, config.nx[pair.alpha]);
	//mapping
	config.xmap(pair.alpha, gr);
	
	cout << "Old grid:" << endl;
	gr.print();
	cout << endl;
}

void BoundGrCalculator::makePotentials() {
	double xConv;
	xConv = 1.0;
	scaleCoo(xConv);
	xConv = 1.0 / xConv;
	vc = make_shared<ScaledPot<Complex>>(\
		pair.vc, eConv, xConv);
	vs = make_shared<ScaledPot<Complex>>(\
		pair.vs, eConv, xConv);
}

void BoundGrCalculator::makeChi() {
	//CPC 182 (2011) 2099–2106 with k = 2
	vector<shared_ptr<Function<Complex, 1>>> phims;
	vector<Int> ls;
	vector<double> epsm;
	Int nob;
	Int l = 0;
	do {
		RadEigenProb tbep(pair, l);
		nob = tbep.numberOfBound();
		for (Int k = 0; k < tbep.numberOfBound(); k++) {
			ls.push_back(l);
			epsm.push_back(tbep.eval[k].real());
			phims.push_back(\
				dynamic_pointer_cast<Function<Complex, 1>>(tbep.getEvec(k)));
		}
		l++;
	} while (nob > 0);
	assert(!phims.empty());

	//regular grid
	double xmin = gr.getLeftmostPoint();
	double xmax = gr.getRightmostPoint();
	Int nPoi = CHI_NPOI_GRID_CALC;
	double h = (xmax - xmin) / (nPoi - 1);
	array<double, 1> xs;
	vector<vector<Complex>> phim_ders;
	phim_ders.resize(phims.size());
	string filename;
	vector<double> f, fd, fd2, fd3, g, gd, gd2;
	vector<Complex> phim, phimd;
	f.resize(nPoi-1); fd.resize(nPoi-1); fd2.resize(nPoi-1);
	g.resize(nPoi-1); gd.resize(nPoi-1); gd2.resize(nPoi-1);
	phim.resize(nPoi - 1); phimd.resize(nPoi - 1); fd3.resize(nPoi-1);
	for (Int k = 0; k < phims.size(); k++) {
		for (Int i = 1; i < nPoi; i++) {
			xs[0] = xmin + h * i;
			fd2[i-1] = fd[i-1] = f[i-1] = \
				(*vs)(xs[0]).real() + (*vc)(xs[0]).real() + \
				ls[k] * (ls[k] + 1) / pow(xs[0], 2.0) - epsm[k];
			phim[i - 1] = (*phims[k])(xs);
			phimd[i - 1] = phims[k]->df(xs);
		}
		firstDeriv(fd, h);
		secondDeriv(fd2, h);
		for (Int i = 1; i < nPoi; i++) {
			xs[0] = xmin + h * i;
			gd2[i - 1] = gd[i - 1] = g[i - 1] = \
				fd2[i - 1] + pow(f[i-1], 2.0);
			fd3[i - 1] = fd2[i - 1];
		}
		firstDeriv(fd3, h);
		firstDeriv(gd, h);
		secondDeriv(gd2, h);

		phim_ders[k].resize(nPoi);
		for (Int i = 1; i < nPoi; i++) {
			xs[0] = xmin + h * i;
			phim_ders[k][i] = \
				( gd2[i-1]+g[i-1]*f[i-1]+ \
						4.0*fd2[i-1]*f[i-1]+2.0*pow(fd[i-1], 2.0))*phim[i-1] + \
				2.0*(gd[i-1]+fd3[i-1]+f[i-1]*fd[i-1])*phimd[i-1];
		}
		phim_ders[k][0] = phim_ders[k][1]; //by continuity
	}
	vector<double> epsx; epsx.resize(nPoi);
	filename = "EPSX_" + config.system + std::to_string(pair.alpha) + "_" + \
		std::to_string((Int)config.xmax[pair.alpha]) + ".dat";
	//ofstream fout(filename);
	for (Int i = 0; i < nPoi; i++) {
		xs[0] = xmin + h * i;
		epsx[i] = 0.0;
		for (Int k = 0; k < phim_ders.size(); k++)
			epsx[i] += pow(abs(phim_ders[k][i]), 2.0);
		epsx[i] = pow(epsx[i], 1.0 / 13.0);
		//fout << xs[0] << "  " << epsx[i] << endl;
	}
	//fout.close();

	//integrate to obtain chi
	vector<double> xi, chim1i;
	xi.reserve(nPoi / 2 + 1); chim1i.reserve(nPoi / 2 + 1);
	xi.push_back(0.0); chim1i.push_back(0.0);
	std::vector<double> vals; vals.resize(3);
	for (Int i = 2; i < nPoi; i += 2) {
		vals[0] = epsx[i - 2]; vals[1] = epsx[i - 1];
		vals[2] = epsx[i];
		xi.push_back((double)i / (nPoi - 1));
		chim1i.push_back(chim1i.back() + simpson(vals, h));
	}
	xi.back() = 1.0;
	for (Int i = 0; i < chim1i.size(); i++)
		chim1i[i] /= chim1i.back();
	filename = "CHIX_" + config.system + std::to_string(pair.alpha) + "_" + \
		std::to_string((Int)config.xmax[pair.alpha]) + ".dat";
	ofstream fout_(filename);
	fout_ << chim1i.size() << endl;
	for (Int i = 0; i < chim1i.size(); i++)
		fout_ << chim1i[i] << "  " << xi[i] << endl;
	fout_.close();

	mapChi(gr, chim1i, xi);
	cout << "New grid:" << endl;
	gr.print();
	cout << endl;
}

void BoundGrCalculator::scaleCoo(double &x) const {
	//physical to reduced coordinates
	pair.reduceCoo(x);
	//scale coordinates due to energy conversion factor
	x /= sqrt(eConv);
}

//BoundGrCalculator::~BoundGrCalculator() { }

/*
shared_ptr<Function<Complex, 1>> RadEigenProb::getEvecDeriv6(const Int i) {
	shared_ptr<AFunction<Complex>> res_ = evec[i]->clone();
	shared_ptr<Function<Complex, 1>> res = \
		dynamic_pointer_cast<Function<Complex, 1>>(res_);

	GenMatrix<Complex> s(nx), d2(nx);
	NumbersList limits;
	NumbersList::const_iterator c;
	double xi, fi; Int j;

	s.fill(zzero); d2.fill(zzero);
	for (Int i = 0; i < nx; i++) {
		xi = cGrid[i];
		basis->getNonzero(xi, limits);
		for (c = limits.begin(); c != limits.end(); c++) {
			j = *c;
			fi = basis->f(xi, j);
			s[i][j] = fi;
			d2[i][j] = basis->dd(xi, j);
		}
	}

	res->coef = d2 * res->coef; s.solve(res->coef);
	res->coef = d2 * res->coef; s.solve(res->coef);
	res->coef = d2 * res->coef; s.solve(res->coef);

	return res;
}
*/