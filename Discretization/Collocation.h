#pragma once

#include "AProjDiscr.h"

template <class T, Int d, class TBas>
class Collocation :
	public AProjDiscr<T, d, TBas> {
public:
	std::array<Grid, d> cGrids;
	Collocation(array<shared_ptr<const ABasis<TBas>>, d> &bases);
	void setCollocGrid(const Int dim, const Grid &grid);
	~Collocation(void) = default;
};

template<class T, Int d, class TBas>
Collocation<T, d, TBas>::Collocation(array<shared_ptr<const ABasis<TBas>>, d> &bases) \
						: AProjDiscr<T, d, TBas>(bases) {
	for (Int k = 0; k < this->bases.size(); k++)	
		cGrids[k] = this->bases[k]->collocGrid();
}

template <class T, Int d, class TBas>
void Collocation<T, d, TBas>::setCollocGrid(const Int dim, const Grid &grid) {
	cGrids[dim] = grid;
}

