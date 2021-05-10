#include "Grid.h"


Grid::Grid(void) {};

Grid::Grid(const vector<double> &points) : poi(points) {
	assert(poi.size() >= 2);
	assert(std::is_sorted(poi.begin(), poi.end()));
}

Grid::Grid(const double lPoint, const double rPoint, const Int nPoints) {
	assert(lPoint <= rPoint);
	assert(nPoints >= 2);
	double h = (rPoint - lPoint) / (nPoints - 1);
	poi.reserve(nPoints);
	for (Int i = 0; i < nPoints; i++) poi.push_back(lPoint + i*h);
}

void Grid::print() const {
	//cout << "Grid   : " << name << endl;
	cout << "Points : " << poi.size() << endl;
	cout << "Bounds : " << double52<double> << poi[0] << \
		" -- " << double52<double> << poi.back() << endl;
	cout << "N Point " << endl;
    for(Int i = 0; i < poi.size(); i++) {
		//cout << i << " " << double135<double> << poi[i] << endl;
		cout << double52<double> << poi[i] << "  ";
    }
	cout << endl;
  }

void Grid::setName(const string &name) {
	this->name = name;
}

string Grid::getName() const {
	return name;
}
