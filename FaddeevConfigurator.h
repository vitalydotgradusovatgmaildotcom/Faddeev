#pragma once

#include "cph.h"
#include "Grid.h"
#include "BinChannel.h"
#include "System3Body.h"

enum FaddeevProblem { \
	bound, scatter, bound_rad, scatter_rad, bound_gr};

class FaddeevConfigurator {
public:
	static const FaddeevConfigurator & getInstance() {
		static FaddeevConfigurator conf;
		return conf;
	}
	System3Body *s3b;
	string system;
	FaddeevProblem problem;
	std::array<double, 3> m, q = {0.0, 0.0, 0.0};
	string idstr;
	Int J, Mproj, tau;
	std::vector<double> energies;
	double eConv = 1.0;
	Int nComp;
	std::array<Int, 3> nx, ny, nz;
	std::array<double, 3> xmax, ymax;
	std::array<string, 3> xGrName;
	std::array<string, 3> yGrName;
	std::vector<Grid> grids;
	std::array<double, 3> x0 = {Inf, Inf, Inf}, \
		y0 = {Inf, Inf, Inf}, nu = {2.01, 2.01, 2.01};
	Int nEig = 1;
	BinChannel chan0;
	Int numPair;
	FaddeevConfigurator(const FaddeevConfigurator &conf) = delete;
	FaddeevConfigurator(FaddeevConfigurator &&conf) = delete;
    FaddeevConfigurator & operator=(const FaddeevConfigurator &conf) = delete;
	FaddeevConfigurator & operator=(FaddeevConfigurator &&conf) = delete;
	void xmap(const Int alpha, Grid &x) const;
	void ymap(const Int alpha, Grid &y) const;
	void zmap(const Int alpha, Grid &z) const;
	~FaddeevConfigurator(void);
protected:
	FaddeevConfigurator(void);
	string getInpFilename();
	void readInpFile(const string &filename);
};

