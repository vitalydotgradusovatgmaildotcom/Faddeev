#include "FaddeevConfigurator.h"
#include "Mappings.h"
#include "SpecialFunctions.h"

FaddeevConfigurator::FaddeevConfigurator(void) {
	readInpFile(getInpFilename());
	//make system
	s3b = new System3Body(system, m, q, idstr, J, Mproj, tau);
}

string FaddeevConfigurator::getInpFilename() {
	std::ifstream f("inputfilename.dat");
	assert(f.is_open());

	string filename;
	getline(f, filename);
	filename = trim(filename);
	filename += ".dat";

	f.close();
	return filename;
}

void FaddeevConfigurator::readInpFile(const string &filename) {
	

	std::ifstream f(filename);
	assert(f.is_open());

	string str;
	while (!f.eof()) {
		getline(f, str);
		if ( str.find(".system") != string::npos ) {
			getline(f, system);
			system = trim(system);
		} else if (str.find(".problem") != string::npos) {
			string problem_;
			getline(f, problem_);
			problem_ = trim(problem_);
			if (problem_ == "bound")
				problem = bound;
			else if (problem_ == "scattering")
				problem = scatter;
			else if (problem_ == "bound radial") {
				problem = bound_rad;
				f >> numPair;
			}
			else if (problem_ == "scatter radial") {
				problem = scatter_rad;
				f >> numPair;
			}
			else if (problem_ == "bound grid") {
				problem = bound_gr;
				f >> numPair;
			}
			else
				assert(false);
		}
		else if ( str.find(".masses") != string::npos )
			f >> m[0] >> m[1] >> m[2];
		else if ( str.find(".charges") != string::npos )
			f >> q[0] >> q[1] >> q[2];
		else if (str.find(".identical") != string::npos) {
			getline(f, idstr);
			idstr = trim(idstr);
		} else if (str.find(".quantum numbers") != string::npos) {
			Int tauPhys; //Physical parity
			f >> J >> Mproj >> tauPhys;
			tau = tauPhys*pow(-1.0, J); //NB!
		} else if ( str.find(".energies") != string::npos ) {
			vector<string> strs;
			getline(f, str);
			split(str, strs);
			if (strs.size() == 1)
				while (str.find_first_not_of(" \t") != string::npos) {
					energies.push_back(atof(str.c_str()));
					getline(f, str);
				}
			else if (strs.size() == 3) {
				energies.push_back(stod(strs[0]));
				double e2 = stod(strs[1]);
				double de = stod(strs[2]);
				double e = energies.back() + de;
				while (e <= e2) {
					energies.push_back(e);
					e += de;
				}
			}
			else
				assert(false);
			assert(energies.size() > 0);
		} else if (str.find(".econv") != string::npos) {
			f >> eConv;
		} else if (str.find(".discretization") != string::npos) {
			//f >> nComp;
			string discr; vector<string> discrs;
			getline(f, discr);
			nComp = atoi(discr.c_str());
			for (Int alpha = 0; alpha < nComp; alpha++) {
				getline(f, discr); split(discr, discrs);
				nx[alpha] = atoi(discrs[0].c_str());
				xmax[alpha] = atof(discrs[1].c_str());
				if (discrs.size() > 2)
					xGrName[alpha] = discrs[2];
				discrs.resize(0);
				getline(f, discr); split(discr, discrs);
				ny[alpha] = atoi(discrs[0].c_str());
				ymax[alpha] = atof(discrs[1].c_str());
				if (discrs.size() > 2)
					yGrName[alpha] = discrs[2];
				discrs.resize(0);
				getline(f, discr);
				nz[alpha] = atoi(discr.c_str());
			}
		}
		else if ( str.find(".cutoff") != string::npos ) {
			for (Int alpha = 0; alpha < 3; alpha++)
				f >> nu[alpha] >> x0[alpha] >> y0[alpha];
		} else if (str.find(".number of eigenvalues") != string::npos) {
			f >> nEig;
		}
		/*
		else if (str.find(".grids") != string::npos) {
			grids.reserve(5);
			string str_, gname;
			Int n; vector<double> points;
			getline(f, str_);
			while (!str_.empty()) {
				gname = str_;
				getline(f, str_);
				n = atoi(str_.c_str());
				points.resize(n);
				for (Int i = 0; i < n; i++) {
					getline(f, str_);
					points[i] = atof(str_.c_str());
				}
				grids.push_back(Grid(points));
				grids.back().setName(gname);
				getline(f, str_);
			}
		}*/
		else if (str.find(".init state") != string::npos) {
			f >> chan0.alpha >> chan0.n >> chan0.l >> chan0.m;
		}
	}

	f.close();
}

void FaddeevConfigurator::xmap(const Int alpha, Grid &x) const {
	if (xGrName[alpha].empty()) {
		//KvitsHu_log(x);
		KvitsHu_pow(x, 4.0); //TODO better KvitsHu_pow_my
	} else { //given grid

		ifstream f;
		f.open(xGrName[alpha]);
		assert(f.is_open());
		Int size;
		std::vector<double> xs, chi;
		f >> size;
		xs.resize(size); chi.resize(size);
		for (Int i = 0; i < size; i++)
			f >> xs[i] >> chi[i];
		f.close();

		mapChi(x, xs, chi);
	}
}

void FaddeevConfigurator::ymap(const Int alpha, Grid &y) const {
	if (yGrName[alpha].empty()) {

		if (system == "HeMol_TTY" || system == "He2Li_KTTY" || \
			system == "He2Na_KTTY" || system == "He2K_KTTY" || \
			system == "He2Rb_KTTY" || system == "He2Cs_KTTY") {
			Roudnev_pow(y, 3.0);
			return;
		}

		//KvitsHu_log(y);
		//KvitsHu_pow(y, 4.0); //TODO better KvitsHu_pow_my

		//mapping for scattering
		//TODO not the best choice I guess
		//mapping must depend on the radius of solution exp decreasing
		Roudnev_pow(y, 1.25);
	} else { //given grid
		assert(false);
		//TODO as in xmap
	}
}

void FaddeevConfigurator::zmap(const Int alpha, Grid &z) const {
	
	Kornev_sin(z);
	//Roudnev_cubic(z, 2.3);
	//Roudnev_tanh(z, 1.25);

//==============================
	/*
	const Int nPoi = 10001; //NB! must be odd for Simpson
	const Int k = 5; //degree of spline
	const Int deg = 3 * z.getNIntervals(); //maximum degree (=num of colloc points)
	//NB! 3 = number of splines per node

	Int Mmin = (1 - tau) / 2, Mmax = J;
	Int nM = Mmax - Mmin + 1;
	std::vector<std::vector<double>> eps_; eps_.resize(nM);
	double h = 1.0 / (nPoi - 1);
	double zpoi, deriv;
	for (Int im = Mmin; im <= Mmax; im++) {
		Int imkp1 = im + k + 1;
		eps_[im-Mmin].resize(nPoi);
		double mult = 1.0, mult2 = 1.0, mult3;
		for (Int imbar = 1; imbar <= imkp1; imbar++)
			mult2 *= (double)(2 * imbar - 1) / (2 * imbar);
		mult2 = sqrt(mult2);
		mult2 *= pow(-1.0, imkp1);
		std::fill(eps_[im-Mmin].begin(), eps_[im - Mmin].end(), 0.0);
		for (Int il = imkp1; il <= im + deg; il++) {
			mult3 = pow(-1.0, k + 1);
			for (Int ik = 1; ik <= k + 1; ik++)
				mult3 *= sqrt(il*(il + 1.0) - (-im - ik) * (-im - ik + 1.0));
			for (Int i = 0; i < nPoi - 1; i++) {
				zpoi = i * h;
				deriv = mult3 * sqrt(2 * PI)*pLegendre(il, imkp1, zpoi) / \
					pow((1 - zpoi)*(1 + zpoi), 0.5*imkp1);
				eps_[im-Mmin][i] += deriv * deriv;
			}
			deriv = mult3 * sqrt(0.5*(2 * il + 1)*mult)*mult2;
			eps_[im-Mmin].back() += deriv * deriv;
			mult *= (double)(il + 1 + imkp1) / (il + 1 - imkp1);
		}
	}

	std::vector<double> eps; eps.resize(nPoi);
	std::fill(eps.begin(), eps.end(), 0.0);
	for (Int i = 0; i < nPoi; i++)
		for (Int im = Mmin; im <= Mmax; im++)
			eps[i] += eps_[im-Mmin][i];

	for (Int i = 0; i < nPoi; i++)
		eps[i] = pow(eps[i], 1.0 / (2 * k + 3));

	ofstream f("EPSZ.dat");
	for (Int i = 0; i < nPoi; i++) {
		zpoi = i * h;
		f << zpoi << "  " << eps[i] << endl;
	}
	f.close();

	std::vector<double> zs;
	std::vector<double> chim1;
	zs.reserve(nPoi / 2 + 1); chim1.reserve(nPoi / 2 + 1);
	zs.push_back(0.0); chim1.push_back(0.0);
	std::vector<double> vals; vals.resize(3);
	for (Int i = 2; i < nPoi; i += 2) {
		vals[0] = eps[i - 2]; vals[1] = eps[i - 1];
		vals[2] = eps[i];
		zs.push_back((double)i / (nPoi - 1));
		chim1.push_back(chim1.back() + simpson(vals, h));
	}
	zs.back() = 1.0;
	for (Int i = 0; i < chim1.size(); i++)
		chim1[i] /= chim1.back();

	f.open("CHIZ.dat");
	for (Int i = 0; i < chim1.size(); i++)
		f << chim1[i] << "  " << zs[i] << endl;
	f.close();

	std::vector<double> chim1_, zs_;
	chim1_.reserve(chim1.size() - 1);
	zs_.reserve(chim1.size() - 1);
	for (Int i = chim1.size() - 1; i > 0; i--) {
		chim1_.push_back(0.5*(-chim1[i] + 1.0));
		zs_.push_back(-zs[i]);
	}
	chim1_.push_back(0.5*(chim1[0] + 1.0));
	zs_.push_back(zs[0]);
	for (Int i = 1; i < chim1.size(); i++) {
		chim1_.push_back(0.5*(chim1[i] + 1.0));
		zs_.push_back(zs[i]);
	}

	//f.open("CHIZ.dat");
	//for (Int i = 0; i < chim1_.size(); i++)
	//	f << chim1_[i] << "  " << zs_[i] << endl;
	//f.close();

	Grid gr(0.0, 1.0, z.getNPoints());
	mapChi(gr, chim1_, zs_);
	for (Int i = 0; i < z.getNPoints(); i++)
		z[i] = gr[i];

	*/
}

FaddeevConfigurator::~FaddeevConfigurator() {
	if (s3b != nullptr)
		delete s3b;
}
