#pragma once

#include "cph.h"

class Grid {
public:
	Grid(void);
	Grid(const vector<double> &points);
	Grid(const double lPoint, const double rPoint, const Int nPoints);
	Grid(const Grid &Grid) = default;
	Grid(Grid &&Grid) = default;
	Grid & operator=(const Grid &rhs) = default;
	Grid & operator=(Grid &&rhs) = default;
	inline Int getNPoints() const;
	inline Int getNIntervals() const; //number of Interval enumerated from 0
	inline Int getInterval(const double x) const;
	inline double & operator[](const Int i);
	inline const double & operator[](const Int i) const;
	inline double getLeftmostPoint() const;
	inline double getRightmostPoint() const;
	void print() const;
	void setName(const string &name);
	string getName() const;
	~Grid(void) = default;
protected:
	string name;
	std::vector<double> poi;	//points
};

inline Int Grid::getNPoints() const {
	return poi.size();
}

inline Int Grid::getNIntervals() const {
	return poi.size() - 1;
}

inline Int Grid::getInterval(const double x) const {
	assert(x >= getLeftmostPoint() && x <= getRightmostPoint());
	if (x == getRightmostPoint())
		return getNIntervals() - 1;
	static vector<double>::const_iterator up;
	up = std::upper_bound(poi.begin(), poi.end(), x);
	return --up - poi.begin();
}

inline double & Grid::operator[](const Int i) {
	return poi[i];
}

inline const double & Grid::operator[](const Int i) const {
	return poi[i];
}

inline double Grid::getLeftmostPoint() const {
	return poi[0];
}

inline double Grid::getRightmostPoint() const {
	return poi.back();
}

