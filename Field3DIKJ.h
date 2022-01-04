#pragma once
#include "Field3D.h"

class Field3DIKJ : public Field3D
{
public:
	Field3DIKJ(int NX, int NY, int NZ) : Field3D(NX, NY, NZ) {};
	Field3DIKJ(int NX, int NY, int NZ, std::vector<double>&& field0) : Field3D(NX, NY, NZ, std::move(field0)) {};
	Field3DIKJ(int NX, int NY, int NZ, const std::vector<double>& field0) : Field3D(NX, NY, NZ, field0) {};
	Field3DIKJ() {}
	~Field3DIKJ() {};

	inline virtual int GetIndex(int i, int j, int k) const override { return i + (k*cNX) + (j*cNX*cNZ); }

	inline double& operator()(int i, int j, int k) override
	{
		return cField[i + (k*cNX) + (j*cNX*cNZ)];
	}

	inline double operator()(int i, int j, int k) const override
	{
		return cField[i + (k*cNX) + (j*cNX*cNZ)];
	}

};