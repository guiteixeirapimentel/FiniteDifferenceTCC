#pragma once
#include <vector>

class Field3D
{
	friend class LSSolver; // para permitir o uso do operador [];
public:
	Field3D(int NX, int NY, int NZ);
	Field3D(int NX, int NY, int NZ, std::vector<double>&& field0);
	Field3D(int NX, int NY, int NZ, const std::vector<double>& field0);
	Field3D();
	~Field3D();
	
	inline virtual double& operator()(int i, int j, int k) = 0;
	inline virtual double operator()(int i, int j, int k) const = 0;

	inline virtual int GetIndex(int i, int j, int k) const = 0;
	 
	inline void clear()
	{
		cField.clear();
		cNX = 0;
		cNY = 0;
		cNZ = 0;
	}
	inline void resize(int nn, int NX, int NY, int NZ)
	{
		cField.resize(nn);
		cNX = NX;
		cNY = NY;
		cNZ = NZ;
	}
	inline int size() const
	{
		return cField.size();
	}
	 
	inline void resize(int nn, int NX, int NY, int NZ, double val)
	{
		cField.resize(nn, val);
		cNX = NX;
		cNY = NY;
		cNZ = NZ;
	}
	 
	inline double* data() { return cField.data(); }
	inline const double* const data() const { return cField.data(); }

protected:
	inline double& operator[](int index) { return cField[index]; }
	inline double operator[](int index) const { return cField[index]; }

public:
	int cNX;
	int cNY;
	int cNZ;

	std::vector<double> cField;
};