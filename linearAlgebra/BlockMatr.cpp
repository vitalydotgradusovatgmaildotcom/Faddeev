#include "BlockMatr.h"

template class BlockMatr<double>;
template class BlockMatr<Complex>;

template <class T>
BlockMatr<T>::BlockMatr( \
	const vector<vector<shared_ptr<AMatrix<T>>>> &blocks) \
	: blocks(blocks) {
	assert(!blocks.empty());
	for (Int ib = 0; ib < blocks.size(); ib++)
		assert(std::find( \
			blocks[ib].begin(), blocks[ib].end(), nullptr) == blocks[ib].end());


	nBlock = blocks.size();
	ns.reserve(nBlock);
	n = 0;
	for (Int ib = 0; ib < nBlock; ib++) {
		assert(blocks[ib].size() == nBlock);
		ns.push_back(blocks[ib][ib]->nrows());
		n += ns.back();
		for (Int jb = 0; jb < nBlock; jb++) {
			assert(blocks[ib][jb]->nrows() == blocks[ib][ib]->nrows());
			assert(blocks[jb][ib]->ncols() == blocks[ib][ib]->ncols());
		}
	}
	m = n;
	
	inverse.resize(nBlock);
	for (Int ib = 0; ib < nBlock; ib++) {
		inverse[ib].resize(nBlock);
		std::fill(inverse[ib].begin(), inverse[ib].end(), false);
	}
}

template <class T>
BlockMatr<T>::BlockMatr(\
	const vector<vector<shared_ptr<AMatrix<T>>>> &blocks, \
		const vector<vector<bool>> inverse) \
	: blocks(blocks), inverse(inverse) {
	assert(!blocks.empty());
	for (Int ib = 0; ib < blocks.size(); ib++)
		assert(std::find(\
			blocks[ib].begin(), blocks[ib].end(), nullptr) == blocks[ib].end());
	assert(blocks.size() == inverse.size());
	for (Int ib = 0; ib < blocks.size(); ib++)
		assert(blocks[ib].size() == inverse[ib].size());

	nBlock = blocks.size();
	ns.reserve(nBlock);
	n = 0;
	for (Int ib = 0; ib < nBlock; ib++) {
		assert(blocks[ib].size() == nBlock);
		ns.push_back(blocks[ib][ib]->nrows());
		n += ns.back();
		for (Int jb = 0; jb < nBlock; jb++) {
			assert(blocks[ib][jb]->nrows() == blocks[ib][ib]->nrows());
			assert(blocks[jb][ib]->ncols() == blocks[ib][ib]->ncols());
			assert(!(inverse[ib][jb] && \
				(blocks[ib][jb]->nrows() != blocks[ib][jb]->ncols())));
		}
	}
	m = n;
}

template <class T>
BlockMatr<T> & BlockMatr<T>::operator=( \
								const BlockMatr<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		blocks = rhs.blocks;
		inverse = rhs.inverse;
		nBlock = rhs.nBlock;
		ns = rhs.ns;
	}
	return *this;
}

template <class T>
BlockMatr<T> & BlockMatr<T>::operator=(const AMatrix<T> &rhs) {
	const BlockMatr<T> &rhs_ = \
		dynamic_cast<const BlockMatr<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
BlockMatr<T> & BlockMatr<T>::operator=(BlockMatr<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	blocks = std::move(rhs.blocks);
	inverse = std::move(rhs.inverse);
	nBlock = rhs.nBlock; rhs.nBlock = 0;
	ns = std::move(rhs.ns);
	return *this;
}

template <class T>
BlockMatr<T> & BlockMatr<T>::operator=(AMatrix<T> &&rhs) {
	BlockMatr<T> && rhs_ = dynamic_cast<BlockMatr<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void BlockMatr<T>::fill(const T val) {
	for (Int ib = 0; ib <= nBlock; ib++)
		for (Int jb = 0; jb <= nBlock; jb++)
			if (!inverse[ib][jb])
				blocks[ib][jb]->fill(val);
			else
				blocks[ib][jb]->fill(1.0 / val);
}

template <class T>
void BlockMatr<T>::resize(const Int nnew) {
	assert(false);
}

template <class T>
void BlockMatr<T>::solve(Vector<T> & rhssol) {
	assert(false);
}

template <class T>
void BlockMatr<T>::print() const {
	cout << "Block matrix with blocks:" << endl;
	for (Int ib = 0; ib < nBlock; ib++) {
		for (Int jb = 0; jb < nBlock; jb++) {
			if (inverse[ib][jb])
				cout << "Block (" << ib << ", " << jb << \
					"): inverse of matrix" << endl;
			else
				cout << "Block (" << ib << ", " << jb << \
					"): matrix" << endl;
			blocks[ib][jb]->print();
		}
	}
}

template <class T>
void BlockMatr<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Block matrix with blocks:" << endl;
	for (Int ib = 0; ib < nBlock; ib++) {
		for (Int jb = 0; jb < nBlock; jb++) {
			if (inverse[ib][jb])
				f << "Block (" << ib << ", " << jb << \
					"): inverse of matrix" << endl;
			else
				f << "Block (" << ib << ", " << jb << "): matrix" << endl;
			blocks[ib][jb]->write(filename);
		}
	}
	f.close();
}

template <class T>
BlockMatr<T> & BlockMatr<T>::operator*=(const T c) {
	for (Int ib = 0; ib <= nBlock; ib++)
		for (Int jb = 0; jb <= nBlock; jb++)
			if (inverse[ib][jb]) {
				assert(c != T());
				*(blocks[ib][jb]) *= 1.0 / c;
			}
			else
				*(blocks[ib][jb]) *= c;
	return *this;
}

template <class T>
Vector<T> BlockMatr<T>::operator*(const Vector<T> &vec) const {
	Int size = vec.size();
	assert(size == m);
	Vector<T> res(n);
	res.fill(T());

	vector<Int> inds; inds.reserve(nBlock);
	inds.push_back(0);
	for (Int ib = 1; ib < nBlock; ib++)
		inds.push_back(inds.back()+ns[ib-1]);
	vector<Vector<T>> vecs = vec.makePartition(inds);
	vector<Vector<T>> ress = res.makePartition(inds);

	Vector<T> tmp;
	for (Int jb = 0; jb < nBlock; jb++) {
		for (Int ib = 0; ib < nBlock; ib++) {
			if (inverse[ib][jb]) {
				tmp = vecs[jb];
				blocks[ib][jb]->solve(tmp);
			}
			else
				tmp = *blocks[ib][jb] * vecs[jb];
			ress[ib] += tmp;
		}
	}
	return res;
}

template <class T>
double BlockMatr<T>::sizeGb() const {
	assert(false);
}