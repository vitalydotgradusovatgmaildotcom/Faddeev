#include "BlockDiagMatr.h"

template class BlockDiagMatr<double>;
template class BlockDiagMatr<Complex>;

template <class T>
BlockDiagMatr<T>::BlockDiagMatr(const std::vector<std::shared_ptr<AMatrix<T>>> &blocks) \
			: blocks(blocks) {
	assert( std::find(blocks.begin(), blocks.end(), nullptr) == blocks.end() );
	assert( !blocks.empty() );

	n = 0;
	ns.resize(blocks.size());
	for (Int ib = 0; ib < blocks.size(); ib++) {
		ns[ib] = blocks[ib]->nrows();
		n += ns[ib];
	}
	m = n;
	inverse.resize(blocks.size());
	std::fill( inverse.begin(), inverse.end(), false );
}

template <class T>
BlockDiagMatr<T>::BlockDiagMatr(const std::vector<std::shared_ptr<AMatrix<T>>> &blocks, \
			const std::vector<bool> &inverse) \
			: blocks(blocks), inverse(inverse) {
	assert( std::find(blocks.begin(), blocks.end(), nullptr) == blocks.end() );
	assert( !blocks.empty() );
	assert(blocks.size() == inverse.size());

	n = 0;
	ns.resize(blocks.size());
	for (Int ib = 0; ib < blocks.size(); ib++) {
		ns[ib] = blocks[ib]->nrows();
		n += ns[ib];
	}
	m = n;
}

template <class T>
BlockDiagMatr<T> & BlockDiagMatr<T>::operator=(const BlockDiagMatr<T> &rhs) {
	if (this != &rhs) {
		n = rhs.n; m = rhs.m;
		blocks = rhs.blocks;
		inverse = rhs.inverse;
		ns = rhs.ns;
	}
	return *this;
}

template <class T>
BlockDiagMatr<T> & BlockDiagMatr<T>::operator=(const AMatrix<T> &rhs) {
	const BlockDiagMatr<T> &rhs_ = dynamic_cast<const BlockDiagMatr<T> &>(rhs);
	this->operator=(rhs_);
	return *this;
}

template <class T>
BlockDiagMatr<T> & BlockDiagMatr<T>::operator=(BlockDiagMatr<T> &&rhs) {
	n = rhs.n; rhs.n = 0;
	m = rhs.m; rhs.m = 0;
	blocks = std::move(rhs.blocks);
	inverse = std::move(rhs.inverse);
	ns = std::move(rhs.ns);
	return *this;
}

template <class T>
BlockDiagMatr<T> & BlockDiagMatr<T>::operator=(AMatrix<T> &&rhs) {
	BlockDiagMatr<T> && rhs_ = dynamic_cast<BlockDiagMatr<T> &&>(rhs);
	this->operator=(std::move(rhs_));
	return *this;
}

template <class T>
void BlockDiagMatr<T>::fill(const T val) {
	for (Int ib = 0; ib <= blocks.size(); ib++)
		if (!inverse[ib])
			blocks[ib]->fill(val);
		else
			blocks[ib]->fill(1.0 / val);
}

template <class T>
void BlockDiagMatr<T>::resize(const Int nnew) {
	assert(false);
}

template <class T>
void BlockDiagMatr<T>::solve(Vector<T> & rhssol) {
	Int size = rhssol.size();
	assert(size == n);

	vector<Int> inds; inds.reserve(blocks.size());
	inds.push_back(0);
	for (Int ib = 1; ib < blocks.size(); ib++)
		inds.push_back(inds.back()+ns[ib-1]);
	vector<Vector<T>> vecs = rhssol.makePartition(inds);

	for (Int ib = 0; ib < blocks.size(); ib++)
		if (inverse[ib]) {
			Vector<T> tmp(*blocks[ib] * vecs[ib]);
			vecs[ib] = tmp;
		}
		else
			blocks[ib]->solve(vecs[ib]);
}

template <class T>
void BlockDiagMatr<T>::print() const {
	cout << "Block diagonal matrix with blocks:" << endl;
	for (Int ib = 0; ib < blocks.size(); ib++) {
		if (inverse[ib])
			cout << "Block " << ib \
				<< ": inverse of matrix" << endl;
		else
			cout << "Block " << ib << ": matrix" << endl;
		blocks[ib]->print();
	}
}

template <class T>
void BlockDiagMatr<T>::write(const string &filename) const {
	ofstream f;
	f.open(filename, ofstream::ate);
	assert(f);
	f << "Block diagonal matrix with blocks:" << endl;
	for (Int ib = 0; ib < blocks.size(); ib++) {
		if (inverse[ib])
			f << "Block " << ib \
				<< ": inverse of matrix" << endl;
		else
			f << "Block " << ib << ": matrix" << endl;
		blocks[ib]->write(filename);
	}
	f.close();
}

template <class T>
BlockDiagMatr<T> & BlockDiagMatr<T>::operator*=(const T c) {
	for (Int ib = 0; ib <= blocks.size(); ib++)
		if (inverse[ib]) {
			assert(c != T());
			*(blocks[ib]) *= 1.0/c;
		} else
			*(blocks[ib]) *= c;
	return *this;
}

template <class T>
Vector<T> BlockDiagMatr<T>::operator*(const Vector<T> &vec) const {
	Int size = vec.size();
	assert(size == n);
	Vector<T> res(vec);

	vector<Int> inds; inds.reserve(blocks.size());
	inds.push_back(0);
	for (Int ib = 1; ib < blocks.size(); ib++)
		inds.push_back(inds.back()+ns[ib -1]);
	vector<Vector<T>> vecs = res.makePartition(inds);

	for (Int ib = 0; ib < blocks.size(); ib++)
		if (inverse[ib])
			blocks[ib]->solve(vecs[ib]);
		else {
			Vector<T> tmp(*blocks[ib] * vecs[ib]);
			vecs[ib] = tmp;
		}

	return res;
}

template <class T>
double BlockDiagMatr<T>::sizeGb() const {
	assert(false);
}
