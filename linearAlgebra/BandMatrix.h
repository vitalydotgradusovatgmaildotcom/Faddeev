#pragma once

#include "AMatrix.h"

template <class T>
class BandMatrix :
	public AMatrix<T> {
public:
	BandMatrix(void);
	~BandMatrix(void);
};

