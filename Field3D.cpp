#include "Field3D.h"

Field3D::Field3D(int NX, int NY, int NZ) 
	:
cNX(NX), cNY(NY), cNZ(NZ) 
{
	cField.resize(NX*NY*NZ, 0.0);
}

Field3D::Field3D(int NX, int NY, int NZ, const std::vector<double>& field0)
	:
	cNX(NX), cNY(NY), cNZ(NZ),
	cField(field0)
{}

Field3D::Field3D(int NX, int NY, int NZ, std::vector<double>&& field0)
	:
	cNX(NX), cNY(NY), cNZ(NZ),
	cField(std::move(field0))
{}

Field3D::Field3D()
	:
	cNX(0),
	cNY(0),
	cNZ(0),
	cField({})
{}

Field3D::~Field3D() {}