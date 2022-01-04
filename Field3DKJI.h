#pragma once
#include "Field3D.h"

class Field3DKJI : public Field3D
{
public:
	Field3DKJI(int NX, int NY, int NZ) : Field3D(NX, NY, NZ) {};
	Field3DKJI(int NX, int NY, int NZ, std::vector<double>&& field0) : Field3D(NX, NY, NZ, std::move(field0)) {};
	Field3DKJI(int NX, int NY, int NZ, const std::vector<double>& field0) : Field3D(NX, NY, NZ, field0) {};
	Field3DKJI() {}
	~Field3DKJI() {};

	virtual int GetIndex(int i, int j, int k) const override { return k + (j*cNZ) + (i * cNZ * cNY); }

	inline double& operator()(int i, int j, int k) override
	{
		return cField[k + (j*cNZ) + (i * cNZ * cNY)];
	}

	inline double operator()(int i, int j, int k) const override
	{
		return cField[k + (j*cNZ) + (i * cNZ * cNY)];
	}
};