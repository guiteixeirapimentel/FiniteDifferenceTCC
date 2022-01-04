#pragma once
#include "Field3D.h"

class Field3DJIK : public Field3D
{
public:
	Field3DJIK(int NX, int NY, int NZ) : Field3D(NX, NY, NZ) {};
	Field3DJIK(int NX, int NY, int NZ, std::vector<double>&& field0) : Field3D(NX, NY, NZ, std::move(field0)) {};
	Field3DJIK(int NX, int NY, int NZ, const std::vector<double>& field0) : Field3D(NX, NY, NZ, field0) {};
	Field3DJIK() {};
	~Field3DJIK() {};

	virtual int GetIndex(int i, int j, int k) const override { return j + (i * cNY) + (k * cNY * cNX); }

	inline double& operator()(int i, int j, int k) override
	{
		return cField[j + (i * cNY) + (k * cNY * cNX)];
	}

	inline double operator()(int i, int j, int k) const override
	{
		return cField[j + (i * cNY) + (k * cNY * cNX)];
	}
};