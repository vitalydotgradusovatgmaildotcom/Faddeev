#include "CompsDiscr.h"


CompsDiscr::CompsDiscr(const vector<shared_ptr<const ADiscretization<Complex, 3>>> &discrs, \
				const Int nComp, const Int Mmin, const Int Mmax, \
					const vector<Int> &solComps) \
			: VectDiscr<Complex, 3>(discrs), nComp(nComp), \
				Mmin(Mmin), Mmax(Mmax), nM(Mmax-Mmin+1), \
					solComps(solComps) {
	assert(nM*solComps.size() == discrs.size());
	makeMap();
}

vector<Int> CompsDiscr::getSolComps() const {
	return solComps;
}

void CompsDiscr::makeMap() {
	comp2sol.reserve(nComp);
	comp2sol.push_back(solComps[0]);
	Int alpha = 0; Int ind = 1;
	while (++alpha < nComp && ind < solComps.size())
		if (alpha != solComps[ind])
			comp2sol.push_back(comp2sol.back());
		else {
			comp2sol.push_back(comp2sol.back() + 1);
			ind++;
		}
	while (alpha < nComp) {
		comp2sol.push_back(comp2sol.back());
		alpha++;
	}
}
