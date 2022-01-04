#include "FDColocated3DField.h"
#include "LSSolver.h"

#include <string>

FDColocated3DField::FDColocated3DField(double rho, double mu, double h, double L, double H, double W, double ufreestream)
	:
	crho(rho),
	cmu(mu),
	ch(h),
	cL(L),
	cH(H),
	cW(W),
	
	cFirstIteration(true),

	cNXU(int(cL / ch) + 1),
	cNYU(int(cW / ch) + 2),
	cNZU(int(cH / ch) + 2),

	cNXV(int(cL / ch) + 2),
	cNYV(int(cW / ch) + 1),
	cNZV(int(cH / ch) + 2),

	cNXW(int(cL / ch) + 2),
	cNYW(int(cW / ch) + 2),
	cNZW(int(cH / ch) + 1),

	cNXP(int(cL / ch) + 2),
	cNYP(int(cW / ch) + 2),
	cNZP(int(cH / ch) + 2),

	cUfreestream(ufreestream)
{
	cUn.resize(cNXU*cNYU*cNZU, 0.0);
	cVn.resize(cNXV*cNYV*cNZV, 0.0);
	cWn.resize(cNXW*cNYW*cNZW, 0.0);
	cPn.resize(cNXP*cNYP*cNZP, 0.0);

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
	cPnp = cPn;

	cPhi = cPn;

	cdphidx = cPn;
	cdphidy = cPn;
	cdphidz = cPn;

	cd2phidx2 = cPn;
	cd2phidy2 = cPn;
	cd2phidz2 = cPn;

	cdUdtn = cUn;
	cdUdt1 = cUn;
	cdUdt2 = cUn;
	cdUdt3 = cUn;

	cdVdtn = cVn;
	cdVdt1 = cVn;
	cdVdt2 = cVn;
	cdVdt3 = cVn;

	cdWdtn = cWn;
	cdWdt1 = cWn;
	cdWdt2 = cWn;
	cdWdt3 = cWn;
}

FDColocated3DField::~FDColocated3DField()
{}

void FDColocated3DField::SetBCFlatPlate()
{
	// freestream in all

	/*for (int j = 0; j < cNY; j++)
	{
		for (int k = 0; k < cNZ; k++)
		{
			for (int i = 0; i < cNX; i++)
			{
				cUn[i + (k * cNX) + (j * cNX * cNZ)] = cUfreestream * ((k*ch + ch)/cH);
				cVn[i + (k * cNX) + (j * cNX * cNZ)] = 0.0;
				cWn[i + (k * cNX) + (j * cNX * cNZ)] = 0.0;
			}
		}
	}*/

	// no slip on z = 0;

	/*for (int j = 1; j < cNY - 1; j++)
	{
		for (int i = 1; i < cNX - 1; i++)
		{
			cUn[i + (0 * cNX) + (j * cNX * cNZ)] = cUfreestream;
			cVn[i + (0 * cNX) + (j * cNX * cNZ)] = 0.0;
			cWn[i + (0 * cNX) + (j * cNX * cNZ)] = 0.0;
		}
	}*/

	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int i = 1; i < cNXU - 1; i++)
		{
			cUn[i + (0 * cNXU) + (j * cNXU * cNZU)] = 0.0;
		}
	}

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
}

void FDColocated3DField::UpdateBCFlatPlate()
{
	// extrapolate outlet (i = NX - 1)

	/*for (int j = 1; j < cNY - 1; j++)
	{
		for (int k = 1; k < cNZ - 1; k++)
		{
			cUn[cNX - 1 + (k * cNX) + (j * cNX * cNZ)] = cUn[cNX - 2 + (k * cNX) + (j * cNX * cNZ)];
			cVn[cNX - 1 + (k * cNX) + (j * cNX * cNZ)] = cVn[cNX - 2 + (k * cNX) + (j * cNX * cNZ)];
			cWn[cNX - 1 + (k * cNX) + (j * cNX * cNZ)] = cWn[cNX - 2 + (k * cNX) + (j * cNX * cNZ)];
		}
	}*/

	// set bc u mom eqt

	// left right
	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			cUn[cNXU - 1 + (k*cNXU) + (j*cNXU * cNZU)] = 0.0;
			cUn[0 + (k*cNXU) + (j*cNXU * cNZU)] = 0.0;
		}
	}

	// top bottom

	for (int j = 1; j < cNYU -1; j++)
	{
		for (int i = 1; i < cNXU-1; i++)
		{
			// bottom
			cUn[i + ((cNZU - 1) * cNXU) + (j*cNXU*cNZU)] = -cUn[i + ((cNZU - 2) * cNXU) + (j*cNXU*cNZU)];

			// top (vel == ufreestream) => (Uk+1 + Uk)/2 = ufreestream
			cUn[i + (0 * cNXU) + (j*cNXU*cNZU)] = 2 * cUfreestream - cUn[i + (1 * cNXU) + (j*cNXU*cNZU)];
		}
	}

	// front back

	for (int k = 1; k < cNZU - 1; k++)
	{
		for (int i = 1; i < cNXU - 1; i++)
		{
			// front
			cUn[i + (k * cNXU) + (0*cNXU*cNZU)] = -cUn[i + (k * cNXU) + (1*cNXU*cNZU)];

			// back
			cUn[i + (k * cNXU) + ((cNYU -1)*cNXU*cNZU)] = -cUn[i + (k * cNXU) + ((cNYU - 2)*cNXU*cNZU)];
		}
	}

	// end u mom bc.

	//set bc v mom eqt

	// left right
	for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			cVn[cNXV - 1 + (k*cNXV) + (j*cNXV * cNZV)] = -cVn[cNZV - 2 + (k*cNXV) + (j*cNXV * cNZV)];
			cVn[0 + (k*cNXV) + (j*cNXV * cNZV)] = -cVn[1 + (k*cNZV) + (j*cNXV * cNZV)];
		}
	}

	// top bottom

	for (int j = 1; j < cNYV -1; j++)
	{
		for (int i = 1; i < cNXV-1; i++)
		{
			cVn[i + ((cNZV - 1) * cNXV) + (j*cNXV*cNZV)] = -cVn[i + ((cNZV - 2) * cNXV) + (j*cNXV*cNZV)];

			// bottom
			cVn[i + (0 * cNXV) + (j*cNXV*cNZV)] = -cVn[i + (1 * cNXV) + (j*cNXV*cNZV)];
		}
	}

	// front back

	for (int k = 1; k < cNZV - 1; k++)
	{
		for (int i = 1; i < cNXV -1; i++)
		{
			// front
			cVn[i + (k * cNXV) + (0 * cNXV*cNZV)] = 0.0;

			// back
			cVn[i + (k * cNXV) + ((cNYV - 1)*cNXV*cNZV)] = 0.0;
		}
	}

	// end bc v mom eqt

	//set bc w mom eqt

	// left right
	for (int j = 1; j < cNYW-1; j++)
	{
		for (int k = 1; k < cNZW-1; k++)
		{
			cWn[cNXW - 1 + (k*cNXW) + (j*cNXW * cNZW)] = -cWn[cNXW - 2 + (k*cNXW) + (j*cNXW * cNZW)];
			cWn[0 + (k*cNXW) + (j*cNXW * cNZW)] = -cWn[1 + (k*cNXW) + (j*cNXW * cNZW)];
		}
	}

	// top bottom

	for (int j = 1; j < cNYW-1; j++)
	{
		for (int i = 1; i < cNXW-1; i++)
		{
			cWn[i + ((cNZW - 1) * cNXW) + (j*cNXW*cNZW)] = 0.0;

			// bottom
			cWn[i + (0 * cNXW) + (j*cNXW*cNZW)] = 0.0;
		}
	}

	// front back

	for (int k = 1; k < cNZW-1; k++)
	{
		for (int i = 1; i < cNXW-1; i++)
		{
			// front
			cWn[i + (k * cNXW) + (0 * cNXW*cNZW)] = -cWn[i + (k * cNXW) + (1 * cNXW*cNZW)];

			// back
			cWn[i + (k * cNXW) + ((cNYW - 1)*cNXW*cNZW)] = -cWn[i + (k * cNXW) + ((cNYW - 2)*cNXW*cNZW)];
		}
	}

	// end bc W mom eqt

	////set bc p eqt

	//// left right
	//for (int j = 1; j < cNYP -1; j++)
	//{
	//	for (int k = 1; k < cNZP-1; k++)
	//	{
	//		cPn[cNXP - 1 + (k*cNXP) + (j*cNXP * cNZP)] = cPn[cNXP - 1 + (k*cNXP) + (j*cNXP * cNZP)];
	//		cPn[0 + (k*cNXP) + (j*cNXP * cNZP)] = cPn[1 + (k*cNXP) + (j*cNXP * cNZP)];
	//	}
	//}

	//// top bottom

	//for (int j = 1; j < cNYP-1; j++)
	//{
	//	for (int i = 1; i < cNXP-1; i++)
	//	{
	//		cPn[i + ((cNZP - 1) * cNXP) + (j*cNXP*cNZP)] = cPn[i + ((cNZP - 2) * cNXP) + (j*cNXP*cNZP)];

	//		// bottom
	//		cPn[i + (0 * cNXP) + (j*cNXP*cNZP)] = cPn[i + (1 * cNXP) + (j*cNXP*cNZP)];
	//	}
	//}

	//// front back

	//for (int k = 1; k < cNZP-1; k++)
	//{
	//	for (int i = 1; i < cNXP-1; i++)
	//	{
	//		// front
	//		cPn[i + (k * cNXP) + (0 * cNXP*cNZP)] = cPn[i + (k * cNXP) + (1 * cNXP*cNZP)];

	//		// back
	//		cPn[i + (k * cNXP) + ((cNYW - 1)*cNXP*cNZP)] = cPn[i + (k * cNXP) + ((cNYP - 2)*cNXP*cNZP)];
	//	}
	//}

	//// end bc p eqt

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
}

int FDColocated3DField::CreatedphidxSystemTDUNUSED(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for(int j = 0; j < NY; j++)
		for(int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{		
				const int index = i + (k * NX) + (j * NX * NZ);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * phi[i + (k * NX) + (j * NX * NZ)] + 4 * phi[i + 1 + (k * NX) + (j * NX * NZ)] - phi[i + 2 + (k * NX) + (j * NX * NZ)]);
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * phi[i + (k * NX) + (j * NX * NZ)] - 4 * phi[i - 1 + (k * NX) + (j * NX * NZ)] + phi[i - 2 + (k * NX) + (j * NX * NZ)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*phi[i - 1 + (k * NX) + (j * NX * NZ)] - (-8 * 0.025 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (-4 * 0.025 + 1.0)*phi[i + 1 + (k * NX) + (j * NX * NZ)] - (-8 * 0.025 / 3.0 + 1.0 / 6.0)*phi[i + 2 + (k * NX) + (j * NX * NZ)] + (-2 * 0.025 / 3.0*phi[i + 3 + (k * NX) + (j * NX * NZ)]));
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*phi[i + 1 + (k * NX) + (j * NX * NZ)] - (8 * 0.09 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (4 * 0.09 + 1.0)*phi[i - 1 + (k * NX) + (j * NX * NZ)] - (8 * 0.09 / 3.0 + 1.0 / 6.0)*phi[i - 2 + (k * NX) + (j * NX * NZ)] + (2 * 0.09 / 3.0*phi[i - 3 + (k * NX) + (j * NX * NZ)]));
				}
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*phi[i - 2 + (k * NX) + (j * NX * NZ)] - Eimp / 2.0*phi[i - 1 + (k * NX) + (j * NX * NZ)] + Eimp / 2.0*phi[i + 1 + (k * NX) + (j * NX * NZ)] + Fimp / 4.0*phi[i + 2 + (k * NX) + (j * NX * NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidySystemTDUNUSED(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * phi[i + (k * NX) + (j * NX * NZ)] + 4 * phi[i + (k * NX) + ((j+1) * NX * NZ)] - phi[i + (k * NX) + ((j + 2) * NX * NZ)]);
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * phi[i + (k * NX) + (j * NX * NZ)] - 4 * phi[i + (k * NX) + ((j - 1) * NX * NZ)] + phi[i + (k * NX) + ((j - 2) * NX * NZ)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*phi[i + (k * NX) + ((j - 1) * NX * NZ)] - (-8 * 0.025 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (-4 * 0.025 + 1.0)*phi[i + (k * NX) + ((j + 1 )* NX * NZ)] - (-8 * 0.025 / 3.0 + 1.0 / 6.0)*phi[i + (k * NX) + ((j + 2) * NX * NZ)] + (-2 * 0.025 / 3.0*phi[i + (k * NX) + ((j + 3) * NX * NZ)]));
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*phi[i + (k * NX) + ((j + 1) * NX * NZ)] - (8 * 0.09 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (4 * 0.09 + 1.0)*phi[i + (k * NX) + ((j - 1) * NX * NZ)] - (8 * 0.09 / 3.0 + 1.0 / 6.0)*phi[i + (k * NX) + ((j - 2) * NX * NZ)] + (2 * 0.09 / 3.0*phi[i + (k * NX) + ((j - 3) * NX * NZ)]));
				}
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*phi[i + (k * NX) + ((j - 2) * NX * NZ)] - Eimp / 2.0*phi[i + (k * NX) + ((j - 1) * NX * NZ)] + Eimp / 2.0*phi[i + (k * NX) + ((j + 1) * NX * NZ)] + Fimp / 4.0*phi[i + (k * NX) + ((j + 2) * NX * NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidzSystemTDUNUSED(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * phi[i + (k * NX) + (j * NX * NZ)] + 4 * phi[i + ((k + 1) * NX) + ((j) * NX * NZ)] - phi[i + ((k + 2) * NX) + ((j) * NX * NZ)]);
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * phi[i + (k * NX) + (j * NX * NZ)] - 4 * phi[i + ((k - 1) * NX) + ((j) * NX * NZ)] + phi[i + ((k - 2) * NX) + ((j) * NX * NZ)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*phi[i + ((k - 1) * NX) + ((j) * NX * NZ)] - (-8 * 0.025 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (-4 * 0.025 + 1.0)*phi[i + ((k + 1) * NX) + ((j)* NX * NZ)] - (-8 * 0.025 / 3.0 + 1.0 / 6.0)*phi[i + ((k + 2) * NX) + ((j) * NX * NZ)] + (-2 * 0.025 / 3.0*phi[i + ((k + 3) * NX) + ((j) * NX * NZ)]));
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*phi[i + ((k + 1) * NX) + ((j) * NX * NZ)] - (8 * 0.09 / 3 + 0.5)*phi[i + (k * NX) + (j * NX * NZ)]
						+ (4 * 0.09 + 1.0)*phi[i + ((k - 1) * NX) + ((j) * NX * NZ)] - (8 * 0.09 / 3.0 + 1.0 / 6.0)*phi[i + ((k - 2) * NX) + ((j) * NX * NZ)] + (2 * 0.09 / 3.0*phi[i + ((k - 3) * NX) + ((j) * NX * NZ)]));
				}
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*phi[i + ((k - 2) * NX) + ((j) * NX * NZ)] - Eimp / 2.0*phi[i + ((k - 1) * NX) + ((j) * NX * NZ)] + Eimp / 2.0*phi[i + ((k + 1) * NX) + ((j) * NX * NZ)] + Fimp / 4.0*phi[i + ((k + 2) * NX) + ((j) * NX * NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::Created2phidx2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i - 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + 
						phi[i + 1 + (k*NX) + (j*NX*NZ)])
						/(ch*ch);
				}
				else*/ if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] + 
						4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)])/(ch*ch);
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
						4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i - 1 + (k*NX) + (j*NX*NZ)])/(ch*ch);
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b/(4.0*ch*ch))*(phi[i + 2 + (k*NX) + (j*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 2 + (k*NX) + (j*NX*NZ)])
						+(a/(ch*ch))*(phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i -1 + (k*NX) + (j*NX*NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::Created2phidy2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j-1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + (k*NX) + ((j+1)*NX*NZ)])
						/ (ch*ch);
				}
				else */if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1 )*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j + 2 )*NX*NZ)] - phi[i + (k*NX) + ((j + 3 )*NX*NZ)]) / (ch*ch);
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1 )*NX*NZ)] +
						4.0 * phi[i + (k*NX) + ((j - 2 )*NX*NZ)] - phi[i + (k*NX) + ((j - 3 )*NX*NZ)]) / (ch*ch);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1 )*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + (k*NX) + ((j - 1 )*NX*NZ)]) / (ch*ch);
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1 )*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + (k*NX) + ((j - 1 )*NX*NZ)]) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + (k*NX) + ((j + 2 )*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 2 )*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + (k*NX) + ((j + 1 )*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 1 )*NX*NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::Created2phidz2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k-1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + ((k+1)*NX) + ((j)*NX*NZ)])
						/ (ch*ch);
				}
				else*/ if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k + 2 )*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3 )*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
						4.0 * phi[i + ((k - 2 )*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3 )*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
						+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
				}				
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 2 )*NX) + ((j)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::CreateSparsePressureEquation(const field3D& dUdx,
	const field3D& dVdy,
	const field3D& dWdz,
	const double dt,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<double>& val,
	std::vector<double>& rhs)const
{
	const int n2 = cNXP * cNYP * cNZP; // Number of points in the grid.
	
	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 13);

	val.clear();
	val.reserve(n2 * 13);

	rhs.clear();
	rhs.resize(n2);

	for (int j = 0; j < cNYP; ++j)
	{
		for (int k = 0; k < cNZP; k++)
		{
			for (int i = 0; i < cNXP; i++) 
			{
				const int index = i + (k * cNXP) + (j*cNXP*cNZP);
					
				if (i == cNXP/2 && j == cNYP/2 && k == cNZP/2)
				{
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (
					// baixo
					(i == 0 && k == 0) || (i == cNXP - 1 && k == 0)
					|| (j == cNYP - 1 && k == 0) || (k == 0 && j == 0) ||

					// cima
					(i == 0 && k == cNZP - 1) || (i == cNXP - 1 && k == cNZP - 1)
					|| (j == cNYP - 1 && k == cNZP - 1) || (k == cNZP - 1 && j == 0) ||
					
					//laterais
					(i == 0 && j == 0) || (i == cNXP - 1 && j == 0)
					|| (i == 0 && j == cNYP - 1) || (i == cNXP - 1 && j == cNYP - 1) 
					)
				{
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}
				else if (i == cNXP - 1)
				{
					// Outlet => Pressure set to zero
					col.push_back(index);
					val.push_back(1.0);

					/*col.push_back(index - 1);
					val.push_back(-1.0);*/
					
					rhs[index] = 0.0;
				}
				else if (k == 0)
				{
					// surface => dPdn = 0
					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					// jkpi
					/*col.push_back(index + (1 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0);*/

					rhs[index] = 0.0;
				}
				else if (k == cNZP - 1)
				{
					// top => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					// jkpi
					/*col.push_back(index + (-1 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0);*/

					rhs[index] = 0.0;
				}
				else if (i == 0)
				{
					// inlet => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					// jkip
					/*col.push_back(index + 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0);*/

					rhs[index] = 0.0;

				}
				else if (j == cNYP - 1)
				{
					// back => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					// jkpi
					/*col.push_back(index + (0 * cNXP) + (-1 * cNXP * cNZP));
					val.push_back(-1.0);*/

					rhs[index] = 0.0;
				}
				else if (j == 0)
				{
					// front => dPdn = 0

					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(1.0);

					// jkpi
					/*col.push_back(index + (0 * cNXP) + (1 * cNXP * cNZP));
					val.push_back(-1.0);*/

					rhs[index] = 0.0;
				}
				else  /*if (i == 1 || i == 2 ||  i == 3 ||
				i == cNXP - 2 || i == cNXP - 3 || i == cNXP - 4 ||
				j == 1 || j == 2 || j == 3 ||
				j == cNYP - 2 || j == cNYP - 3 || j == cNYP - 4 ||
				k == 1 || k == 2 || k == 3 ||
				k == cNZP - 2 || k == cNZP - 3 || k == cNZP - 4)*/
				{
					// reduz o stencil utilizado para 2
					
					double valijk = 0.0;
					
					// jkmi
					col.push_back(index + (-1 * cNXP) + (0 * cNXP * cNZP));
					if (k == 1)
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkpi
					col.push_back(index + (1 * cNXP) + (0 * cNXP * cNZP));
					if (k == cNZP - 2)
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkim
					col.push_back(index - 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					if (i == 1)
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

										
					// jkip
					col.push_back(index + 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					if (i == cNXP - 2)
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

					// jmki
					col.push_back(index + (0 * cNXP) + (-1 * cNXP * cNZP));
					if (j == 1)
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jpki
					col.push_back(index + (0 * cNXP) + (1 * cNXP * cNZP));
					if (j == cNYP - 2)
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(valijk);


					rhs[index] = ch*ch*(crho / dt)*(
						dUdx[i - 1 + (k*(cNXU-1)) + (j*(cNXU - 1)*cNZU)] +
						dVdy[i + (k*cNXV) + ((j-1)*cNXV *cNZV)] +
						dWdz[i + ((k-1)*cNXW) + (j*cNXW *(cNZW-1))]
						);
				}
				/*else if (i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					// DENTRO DO SOLIDO OU EXATAMENTE NA PAREDE DO SÓLIDO ( SETA PCORR = 0 )
					col.push_back(index);
					val.push_back(1.0);

					rhs[index] = 0.0;
				}*/
				/*else
				{
					// [i + (k * cNX) + (j * cNX * cNZ)] = [index global]
					
					// jkmi
					col.push_back(index + (-1 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkmmi
					col.push_back(index + (-2 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0 / (12.0));
					

					// jkim
					col.push_back(index - 1 + (0 * cNXP) + (0 * cNXP * cNZP));					
					val.push_back(16.0/(12.0));

					// jkimm
					col.push_back(index - 2 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0/(12.0));

					// jkip
					col.push_back(index + 1 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkipp
					col.push_back(index + 2 + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0 / (12.0));
					

					// jkpi
					col.push_back(index + (1 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jkppi
					col.push_back(index + (2 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-1.0/(12.0));

					// jmki
					col.push_back(index + (0 * cNXP) + (-1 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jmmki
					col.push_back(index + (0 * cNXP) + (-2 * cNXP * cNZP));
					val.push_back(-1.0 / (12.0));
					

					// jpki
					col.push_back(index + (0 * cNXP) + (1 * cNXP * cNZP));
					val.push_back(16.0/(12.0));

					// jppki
					col.push_back(index + (0 * cNXP) + (2 * cNXP * cNZP));
					val.push_back(-1.0/(12.0));


					// jki
					col.push_back(index + (0 * cNXP) + (0 * cNXP * cNZP));
					val.push_back(-30.0*3.0 / (12.0));
					
					rhs[index] = ch * ch*(crho / dt)*(
						dUdx[i - 1 + (k*(cNXU - 1)) + (j*(cNXU - 1)*cNZU)] +
						dVdy[i + (k*cNXV) + ((j - 1)*cNXV *cNZV)] +
						dWdz[i + ((k - 1)*cNXW) + (j*cNXW *(cNZW - 1))]
						);
				}*/

				ptr.push_back(col.size());
			}
		}
	}

	return n2;
}

void FDColocated3DField::CalcdphidtUmom(field3D& dphidt, 
	LSSolver& lssolver, 
	std::vector<double>& aa, 
	std::vector<double>& bb, 
	std::vector<double>& cc)
{
	int nn = 0;

	// interpolate v in x
	nn = CreateInterpolateXSystemTD(cVnp, cNXV, cNYV, cNZV, aa, bb, cc, cdphidx);

	field3D vnpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);
	vnpInterpolated = cdphidx;

	// interpolate w in x

	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNYW, cNZW, aa, bb, cc, cdphidx);

	field3D wnpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);
	wnpInterpolated = cdphidx;

	// interpolate u in x, y, z

	nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNYU, cNZU, aa, bb, cc, cphiInterpolateX);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateX, nn, aa, bb, cc);

	nn = CreateInterpolateYSystemTD(cUnp, cNXU, cNYU, cNZU, aa, bb, cc, cphiInterpolateY);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateY, nn, aa, bb, cc);

	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNYU, cNZU, aa, bb, cc, cphiInterpolateZ);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateZ, nn, aa, bb, cc);
	
	// u mom eq

	// calc duudx
	nn = CreatedphidxUConvSystemTD(aa, bb, cc, cdphidx);

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	// calc duvdy
	nn = CreatedphidyUConvSystemTD(vnpInterpolated, aa, bb, cc, cdphidy);

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	// calc duwdz
	nn = CreatedphidzUConvSystemTD(wnpInterpolated, aa, bb, cc, cdphidz);

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cd2phidx2);

	lssolver.SolveTridiagonalDestructive(cd2phidx2, nn, aa, bb, cc);

	// calc d2udy2
	nn = Created2phidy2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cd2phidy2);

	lssolver.SolveTridiagonalDestructive(cd2phidy2, nn, aa, bb, cc);

	// calc d2udz2
	nn = Created2phidz2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, cd2phidz2);

	lssolver.SolveTridiagonalDestructive(cd2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				dphidt[i + (k*cNXU) + (j* cNXU* cNZU)] =
					-1.0*(cdphidx[i + (k* cNXU) + (j* cNXU* cNZU)] +
						 cdphidy[i + (k* cNXU) + (j* cNXU* cNZU)] +
						 cdphidz[i + (k* cNXU) + (j* cNXU* cNZU)])
					+ (cmu / crho)*(cd2phidx2[i + (k* cNXU) + (j* cNXU* cNZU)]
						+  cd2phidy2[i + (k* cNXU) + (j* cNXU* cNZU)] +
						 cd2phidz2[i + (k* cNXU) + (j* cNXU* cNZU)]);
			}
		}
	}
}

void FDColocated3DField::CalcdphidtVmom(field3D& dphidt,
	LSSolver& lssolver,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	// interpolate u in y
	nn = CreateInterpolateYSystemTD(cUnp, cNXU, cNYU, cNZU, aa, bb, cc, cdphidy);

	field3D unpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	unpInterpolated = cdphidy;

	// interpolate w in y

	nn = CreateInterpolateYSystemTD(cWnp, cNXW, cNYW, cNZW, aa, bb, cc, cdphidy);

	field3D wnpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	wnpInterpolated = cdphidy;

	// interpolate v in x, y, z

	nn = CreateInterpolateXSystemTD(cVnp, cNXV, cNYV, cNZV, aa, bb, cc, cphiInterpolateX);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateX, nn, aa, bb, cc);

	nn = CreateInterpolateYSystemTD(cVnp, cNXV, cNYV, cNZV, aa, bb, cc, cphiInterpolateY);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateY, nn, aa, bb, cc);

	nn = CreateInterpolateZSystemTD(cVnp, cNXV, cNYV, cNZV, aa, bb, cc, cphiInterpolateZ);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateZ, nn, aa, bb, cc);

	// v mom eq

	// calc duvdx
	nn = CreatedphidxVConvSystemTD(unpInterpolated, aa, bb, cc, cdphidx);

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	// calc dvvdy
	nn = CreatedphidyVConvSystemTD(aa, bb, cc, cdphidy);

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	// calc dvwdz
	nn = CreatedphidzVConvSystemTD(wnpInterpolated, aa, bb, cc, cdphidz);

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	// calc d2vdx2
	nn = Created2phidx2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cd2phidx2);

	lssolver.SolveTridiagonalDestructive(cd2phidx2, nn, aa, bb, cc);

	// calc d2vdy2
	nn = Created2phidy2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cd2phidy2);

	lssolver.SolveTridiagonalDestructive(cd2phidy2, nn, aa, bb, cc);

	// calc d2vdz2
	nn = Created2phidz2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, cd2phidz2);

	lssolver.SolveTridiagonalDestructive(cd2phidz2, nn, aa, bb, cc);

	// get dvdt

	for (int j = 0; j < cNYV; j++)
	{
		for (int k = 0; k < cNZV; k++)
		{
			for (int i = 0; i < cNXV; i++)
			{
				dphidt[i + (k*cNXV) + (j* cNXV* cNZV)] =
					-1.0*(cdphidx[i + (k* cNXV) + (j* cNXV* cNZV)] +
						cdphidy[i + (k* cNXV) + (j* cNXV* cNZV)] +
						cdphidz[i + (k* cNXV) + (j* cNXV* cNZV)])
					+ (cmu / crho)*(cd2phidx2[i + (k* cNXV) + (j* cNXV* cNZV)]
						+ cd2phidy2[i + (k* cNXV) + (j* cNXV* cNZV)] +
						cd2phidz2[i + (k* cNXV) + (j* cNXV* cNZV)]);
			}
		}
	}
}

void FDColocated3DField::CalcdphidtWmom(field3D& dphidt,
	LSSolver& lssolver,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	// interpolate u in z
	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNYU, cNZU, aa, bb, cc, cdphidz);

	field3D unpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	unpInterpolated = cdphidz;

	// interpolate v in z

	nn = CreateInterpolateZSystemTD(cVnp, cNXV, cNYV, cNZV, aa, bb, cc, cdphidz);

	field3D vnpInterpolated;

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	vnpInterpolated = cdphidz;

	// interpolate w in x, y, z

	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNYW, cNZW, aa, bb, cc, cphiInterpolateX);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateX, nn, aa, bb, cc);

	nn = CreateInterpolateYSystemTD(cWnp, cNXW, cNYW, cNZW, aa, bb, cc, cphiInterpolateY);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateY, nn, aa, bb, cc);

	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, aa, bb, cc, cphiInterpolateZ);

	lssolver.SolveTridiagonalDestructive(cphiInterpolateZ, nn, aa, bb, cc);

	// w mom eq

	// calc duwdx
	nn = CreatedphidxWConvSystemTD(unpInterpolated, aa, bb, cc, cdphidx);

	lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

	// calc dvwdy
	nn = CreatedphidyWConvSystemTD(vnpInterpolated,aa, bb, cc, cdphidy);

	lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

	// calc dwwdz
	nn = CreatedphidzWConvSystemTD(aa, bb, cc, cdphidz);

	lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cd2phidx2);

	lssolver.SolveTridiagonalDestructive(cd2phidx2, nn, aa, bb, cc);

	// calc d2udy2
	nn = Created2phidy2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cd2phidy2);

	lssolver.SolveTridiagonalDestructive(cd2phidy2, nn, aa, bb, cc);

	// calc d2udz2
	nn = Created2phidz2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, cd2phidz2);

	lssolver.SolveTridiagonalDestructive(cd2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 0; j < cNYW; j++)
	{
		for (int k = 0; k < cNZW; k++)
		{
			for (int i = 0; i < cNXW; i++)
			{
				dphidt[i + (k*cNXW) + (j* cNXW* cNZW)] =
					-1.0*(cdphidx[i + (k* cNXW) + (j* cNXW* cNZW)] +
						cdphidy[i + (k* cNXW) + (j* cNXW* cNZW)] +
						cdphidz[i + (k* cNXW) + (j* cNXW* cNZW)])
					+ (cmu / crho)*(cd2phidx2[i + (k* cNXW) + (j* cNXW* cNZW)]
						+ cd2phidy2[i + (k* cNXW) + (j* cNXW* cNZW)] +
						cd2phidz2[i + (k* cNXW) + (j* cNXW* cNZW)]);
			}
		}
	}
}


void FDColocated3DField::SaveIKCutToCSV(const std::string& filename, int j)
{
	std::ofstream arqU(filename + "U.csv");
	
	for (int k = 0; k < cNZU; k++)
	{
		for(int i = 0; i < cNXU; i++)
		{
			arqU << cUnp[i + (k*cNXU) + (j * cNXU*cNZU)] << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "V.csv");

	for (int k = 0; k < cNZV; k++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqV << cVnp[i + (k*cNXV) + (j * cNXV*cNZV)] << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "W.csv");

	for (int k = 0; k < cNZW; k++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqW << cWnp[i + (k*cNXW) + (j * cNXW*cNZW)] << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "P.csv");

	for (int k = 0; k < cNZP; k++)
	{
		for (int i = 0; i < cNXP; i++)
		{
			arqP << cPn[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";
		}

		arqP << "\n";
	}

	arqP.close();

	/*std::ofstream arqDivU(filename + "divU.csv");

	for (int k = 0; k < cNZU; k++)
	{
		for (int i = 0; i < (cNXU-1); i++)
		{
			arqDivU << cdphidx[i + (k*(cNXU - 1)) + (j * (cNXU - 1)*cNZU)]  << "\t";
		}

		arqDivU << "\n";
	}

	arqDivU.close();

	std::ofstream arqDivV(filename + "divV.csv");

	for (int k = 0; k < cNZV; k++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqDivV << cdphidy[i + (k* cNXV) + ((j) * cNXV *cNZV)] << "\t";
		}

		arqDivV << "\n";
	}

	arqDivV.close();

	std::ofstream arqDivW(filename + "divW.csv");

	for (int k = 0; k < cNZW-1; k++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqDivW << cdphidz[i + (k* cNXW) + (j * cNXW *(cNZW-1))] << "\t";
		}

		arqDivW << "\n";
	}

	arqDivW.close();*/

}

void FDColocated3DField::SaveJKCutToCSV(const std::string& filename, int i)
{
	std::ofstream arqU(filename + "TransU.csv");

	for (int j = 0; j < cNYU; j++)
	{
		for (int k = 0; k < cNZU; k++)
		{
			arqU << cUnp[i + (k*cNXU) + (j * cNXU*cNZU)] << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "TransV.csv");

	for (int j = 0; j < cNYV; j++)
	{
		for (int k = 0; k < cNZV; k++)
		{
			arqV << cVnp[i + (k*cNXV) + (j * cNXV*cNZV)] << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "TransW.csv");

	for (int j = 0; j < cNYW; j++)
	{
		for (int k = 0; k < cNZW; k++)
		{
			arqW << cWnp[i + (k*cNXW) + (j * cNXW*cNZW)] << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "TransP.csv");

	for (int j = 0; j < cNYP; j++)
	{
		for (int k = 0; k < cNZP; k++)
		{
			arqP << cPn[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";
		}

		arqP << "\n";
	}

	arqP.close();

	//std::ofstream arqDiv(filename + "Transdiv.csv");

	//for (int j = 0; j < cNYP; j++)
	//{
	//	for (int k = 0; k < cNZP; k++)
	//	{
	//		/*arqDiv << cdphidx[i + (k*cNXP) + (j * cNXP*cNZP)] + cdphidy[i + (k*cNXP) + (j * cNXP*cNZP)] +
	//			cdphidz[i + (k*cNXP) + (j * cNXP*cNZP)] << "\t";*/
	//	}

	//	arqDiv << "\n";
	//}

	//arqDiv.close();
}

int FDColocated3DField::CreatedphidxUConvSystemTD(
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXU * cNYU * cNZU;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYU; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = i + (k * cNXU) + (j * cNXU * cNZU);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					
					rhs[index] = 0.0;
				}
				else if (i == cNXU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == 1 || i == 2 || i == cNXU - 2 || i == cNXU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uipjk = cphiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					const double uimjk = cphiInterpolateX[i - 1 + (k * (cNXU-1)) + (j * (cNXU - 1) * cNZU)];

					rhs[index] = ((uipjk * uipjk)
						- (uimjk*uimjk)
						) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uipjk = 0.5*(
						cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] + cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)]
						);

					const double uimjk = 0.5*(
						cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] + cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)]
						);

					rhs[index] = ((uipjk * uipjk)
						- (uimjk*uimjk)
						) / 
						(ch);
				}*/
				/*else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					
					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const double uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uippjk = cUnp[i + 2 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = 0.5 / ch * (-3 * uijk * uijk +
						4 * uipjk * uipjk -
						 uippjk* uippjk);
				}
				else if (i == cNXU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const double uimjk = cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uimmjk = cUnp[i - 2 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = 0.5 / ch * (3 * uijk*uijk -
						4 * uimjk*uimjk +
						uimmjk*uimmjk);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uimjk = cUnp[i -1 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const double uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uippjk = cUnp[i + 2 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uipppjk = cUnp[i + 3 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = (1.0 / ch) * ((-2 * 0.025 / 3 - 1.0 / 3.0)*uimjk*uimjk -
						(-8 * 0.025 / 3 + 0.5)*uijk*uijk
						+ (-4 * 0.025 + 1.0)*uipjk*uipjk -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*uippjk +
						(-2 * 0.025 / 3.0*uipppjk*uipppjk));
				}
				else if (i == cNXU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
					const double uimjk = cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uimmjk = cUnp[i - 2 + (k * cNXU) + (j * cNXU * cNZU)];
					const double uimmmjk = cUnp[i - 3 + (k * cNXU) + (j * cNXU * cNZU)];

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*uipjk -
						(8 * 0.09 / 3 + 0.5)*uijk
						+ (4 * 0.09 + 1.0)*uimjk -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*uimmjk +
						(2 * 0.09 / 3.0*uimmmjk));
				}*/
				else
				{
					
					const double uipjk = cphiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
					const double uippjk = cphiInterpolateX[i + 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					const double uimjk = cphiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
					const double uimmjk = cphiInterpolateX[i - 2 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (
						-Fimp / 4.0 * uimmjk * uimmjk
						- Eimp / 2.0*uimjk*uimjk 
						+ Eimp / 2.0*uipjk*uipjk 
						+ Fimp / 4.0* uippjk * uippjk);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidyUConvSystemTD(
	const field3D& V,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXU * cNYU * cNZU;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYU; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = i + (k * cNXU) + (j * cNXU * cNZU);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == cNYU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == 1 || j == 2 || j == cNYU - 2 || j == cNYU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijmk = cphiInterpolateY[i + (k *cNXU) + ((j-1) * cNXU * cNZU)];
					const double uijpk = cphiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];

					const double vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
					const double vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];
					
					rhs[index] = ((uijpk * vijpk)
						- (uijmk * vijmk)) / ch;
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijpk = 0.5*(
						cUnp[i + (k * cNXU) + ((j) * cNXU * cNZU)] + cUnp[i + (k * cNXU) + ((j + 1) * cNXU * cNZU)]
						);

					const double uijmk = 0.5*(
						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)]
						);

					rhs[index] = ((uijpk * V[i + (k * (cNXV-1)) + ((j) * (cNXV - 1) * cNZV)])
						- (uijmk * V[i  + (k * (cNXV - 1)) + ((j-1) * (cNXV - 1) * cNZV)])) /ch;
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijk = cUnp[i + (k*cNXU) + (j*cNXU*cNZU)];
					const double uijpk = cUnp[i + (k*cNXU) + ((j+1)*cNXU*cNZU)];
					const double uijppk = cUnp[i + (k*cNXU) + ((j + 2)*cNXU*cNZU)];

					const double vijk = V[i + (k*(cNXV - 1)) + (j * (cNXV - 1) * cNXV)];
					const double vijpk = V[i + (k*(cNXV - 1)) + ((j+1) * (cNXV - 1) * cNXV)];
					const double vijppk = V[i + (k*(cNXV - 1)) + ((j+2) * (cNXV - 1) * cNXV)];

					rhs[index] = 0.5 / ch * (-3 * uijk * V[i + (k * cNXU) + (j * cNXU * cNZV)] +
						4 * uijpk * V[i + (k * cNXU) + ((j + 1) * cNXU * cNZV)] -
						uijppk * V[i + (k * cNXU) + ((j + 2) * cNXU * cNZV)]);
				}
				else if (j == cNYU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + ((j-1) * cNXU * cNZV)] -
						4 * cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 2) * cNXU * cNZV)] +
						cUnp[i + (k * cNXU) + ((j - 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 3) * cNXU * cNZV)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 1) * cNXU * cNZV)] -
						(-8 * 0.025 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + (j * cNXU * cNZV)]
						+ (-4 * 0.025 + 1.0)*cUnp[i + (k * cNXU) + ((j + 1)* cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 1)* cNXU * cNZV)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cUnp[i + (k * cNXU) + ((j + 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 2) * cNXU * cNZV)] +
						(-2 * 0.025 / 3.0*cUnp[i + (k * cNXU) + ((j + 3) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 3) * cNXU * cNZV)]));
				}
				else if (j == cNYU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cUnp[i + (k * cNXU) + ((j + 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j) * cNXU * cNZV)] -
						(8 * 0.09 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + ((j-1) * cNXU * cNZV)]
						+ (4 * 0.09 + 1.0)*cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 2) * cNXU * cNZV)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cUnp[i + (k * cNXU) + ((j - 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 3) * cNXU * cNZV)] +
						(2 * 0.09 / 3.0*cUnp[i + (k * cNXU) + ((j - 3) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 4) * cNXU * cNZV)]));
				}*/
				else
				{
					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					const double uijmk = cphiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
					const double uijmmk = cphiInterpolateY[i + (k *cNXU) + ((j - 2) * cNXU * cNZU)];

					const double uijpk = cphiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];
					const double uijppk = cphiInterpolateY[i + (k*cNXU) + ((j+1)* cNXU * cNZU)];

					const double vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
					const double vijppk = V[i + (k * (cNXV - 1)) + ((j+1) * (cNXV - 1) * cNZV)];

					const double vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];
					const double vijmmk = V[i + (k * (cNXV - 1)) + ((j - 2) * (cNXV - 1) * cNZV)];

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*uijmmk * vijmmk -
						Eimp / 2.0*uijmk * vijmk +
						Eimp / 2.0*uijpk * vijpk +
						Fimp / 4.0*uijppk * vijppk);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidzUConvSystemTD(
	const field3D& W,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXU * cNYU * cNZU;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYU; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = i + (k * cNXU) + (j * cNXU * cNZU);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == cNZU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == 1 || k == 2 || k == cNZU - 2 || k == cNZU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijkp = cphiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU-1))];
					const double uijkm = cphiInterpolateZ[i + ((k-1)*cNXU) + (j*cNXU*(cNZU-1))];

					const double wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const double wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					rhs[index] = ((uijkp * wijkp)
						- (uijkm * wijkm)) /
						(ch);
				}
				/*
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					
					const double uijkp = 0.5*(
						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + ((k+1) * cNXU) + ((j) * cNXU * cNZU)]
						);

					const double uijkm = 0.5*(
						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + ((k-1) * cNXU) + ((j) * cNXU * cNZU)]
						);

					rhs[index] = ((uijkp * W[i + ((k) * (cNXW-1)) + ((j)* (cNXW - 1) * cNZW)])
						- (uijkm * W[i + ((k-1) * (cNXW - 1)) + (j  * (cNXW - 1) * cNZW)])) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)] +
						4 * cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
						cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)]);
				}
				else if (k == cNZU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)] -
						4 * cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] +
						cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)]
						+ (-4 * 0.025 + 1.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)] +
						(-2 * 0.025 / 3.0*cUnp[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZW)]));
				}
				else if (k == cNZU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)]
						+ (4 * 0.09 + 1.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)] +
						(2 * 0.09 / 3.0*cUnp[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 4) * cNXU) + ((j)* cNXU * cNZW)]));
				}*/
				else
				{
					const double uijkp = cphiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
					const double uijkpp = cphiInterpolateZ[i + ((k+1)*cNXU) + (j*cNXU*(cNZU - 1))];

					const double uijkm = cphiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];
					const double uijkmm = cphiInterpolateZ[i + ((k - 2)*cNXU) + (j*cNXU*(cNZU - 1))];

					const double wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const double wijkpp = W[i + ((k+1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					const double wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
					const double wijkmm = W[i + ((k - 2) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (
						- Fimp / 4.0* uijkmm * wijkmm 
						- Eimp / 2.0* uijkm * wijkm 
						+ Eimp / 2.0*uijkp * wijkp 
						+ Fimp / 4.0*uijkpp * wijkpp);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidxVConvSystemTD(
	const field3D& U,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXV * cNYV * cNZV;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYV; j++)
		for (int k = 0; k < cNZV; k++)
			for (int i = 0; i < cNXV; i++)
			{
				const int index = i + (k * cNXV) + (j * cNXV * cNZV);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == cNXV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == 1 || i == 2 || i == cNXV - 2 || i == cNXV - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double vipjk = cphiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const double vimjk = cphiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const double uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
					const double uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];

					rhs[index] = ((vipjk * uipjk)
						- (vimjk * uimjk)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double vipjk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + 1 + ((k) * cNXV) + ((j)* cNXV * cNZV)]
						);

					const double vimjk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i - 1 + ((k) * cNXV) + ((j)* cNXV * cNZV)]
						);

					rhs[index] = ((vipjk * U[i + (k * cNXU) + ((j)* cNXU  * cNZU)])
						- (vimjk * U[i  - 1 + (k * cNXU) + (j  * cNXU * cNZU)])) /
						(ch);
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)] +
						4 * cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZU)] -
						cVnp[i + 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZU)]);
				}
				else if (i == cNXV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1+ (k * cNXU) + (j * cNXU * cNZU)] -
						4 * cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZU)] +
						cVnp[i - 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZU)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZU)] -
						(-8 * 0.025 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)]
						+ (-4 * 0.025 + 1.0)*cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZU)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cVnp[i + 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZU)] +
						(-2 * 0.025 / 3.0*cVnp[i + 3 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZU)]));
				}
				else if (i == cNXV - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)] -
						(8 * 0.09 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZU)]
						+ (4 * 0.09 + 1.0)*cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZU)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cVnp[i - 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZU)] +
						(2 * 0.09 / 3.0*cVnp[i - 3 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZU)]));
				}*/
				else
				{
					const double vipjk = cphiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const double vippjk = cphiInterpolateX[i + 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const double vimjk = cphiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
					const double vimmjk = cphiInterpolateX[i - 2 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];

					const double uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
					const double uippjk = U[i + 1 + (k *cNXU) + (j*cNXU*cNZU)];

					const double uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];
					const double uimmjk = U[i - 2 + (k *cNXU) + (j*cNXU*cNZU)];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*vimmjk * uimmjk -
						Eimp / 2.0*vimjk * uimjk +
						Eimp / 2.0*vipjk * uipjk +
						Fimp / 4.0*vippjk * uippjk);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidyVConvSystemTD(
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXV * cNYV * cNZV;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYV; j++)
		for (int k = 0; k < cNZV; k++)
			for (int i = 0; i < cNXV; i++)
			{
				const int index = i + (k * cNXV) + (j * cNXV * cNZV);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == cNYV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == 1 || j == 2 || j == cNYV - 2 || j == cNYV - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double vijpk = cphiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];

					const double vijmk = cphiInterpolateY[i + (k * cNXV) + ((j-1)*cNXV * cNZV)];

					rhs[index] = ((vijpk*vijpk)
						- (vijmk*vijmk)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double vijpk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k)* cNXV) + ((j+1)* cNXV * cNZV)]
						);

					const double vijmk = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k) * cNXV) + ((j-1)* cNXV * cNZV)]
						);

					rhs[index] = ((vijpk*vijpk)
						- (vijmk*vijmk)) /
						(ch);
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] +
						4 * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] -
						cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)]);
				}
				else if (j == cNYV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] -
						4 * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] +
						cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] -
						(-8 * 0.025 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)]
						+ (-4 * 0.025 + 1.0)*cVnp[i + (k * cNXV) + ((j + 1)* cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1)* cNXV * cNZV)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] +
						(-2 * 0.025 / 3.0*cVnp[i + (k * cNXV) + ((j + 3) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 3) * cNXV * cNZV)]));
				}
				else if (j == cNYV - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] -
						(8 * 0.09 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)]
						+ (4 * 0.09 + 1.0)*cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] +
						(2 * 0.09 / 3.0*cVnp[i + (k * cNXV) + ((j - 3) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 3) * cNXV * cNZV)]));
				}*/
				else
				{
					const double vijpk = cphiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];
					const double vijppk = cphiInterpolateY[i + (k * cNXV) + ((j+1)*cNXV * cNZV)];

					const double vijmk = cphiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];
					const double vijmmk = cphiInterpolateY[i + (k * cNXV) + ((j - 2)*cNXV * cNZV)];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*vijmmk * vijmmk -
						Eimp / 2.0*vijmk * vijmk +
						Eimp / 2.0*vijpk * vijpk +
						Fimp / 4.0*vijppk * vijppk);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidzVConvSystemTD(
	const field3D& W,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXV * cNYV * cNZV;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYV; j++)
		for (int k = 0; k < cNZV; k++)
			for (int i = 0; i < cNXV; i++)
			{
				const int index = i + (k * cNXV) + (j * cNXV * cNZV);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == cNZV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == 1 || k == 2 || k == cNZV - 2 || k == cNZV - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double vijkp = cphiInterpolateZ[i + ((k) * cNXV) + (j* cNXV * (cNZV-1))];

					const double vijkm = cphiInterpolateZ[i + ((k-1) * cNXV) + (j* cNXV * (cNZV - 1))];

					const double wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
					const double wijkm = W[i + ((k-1)* cNXW) + (j*cNXW*cNZW)];

					rhs[index] = ((vijkp * wijkp)
						- (vijkm * wijkm)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double vijkp = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k+1)* cNXV) + ((j)* cNXV * cNZV)]
						);

					const double vijkm = 0.5*(
						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)]
						);

					rhs[index] = ((vijkp * W[i + ((k) * cNXW) + ((j)* cNXW * cNZW)])
						- (vijkm * W[i + ((k - 1) * cNXW) + ((j)  * cNXW * cNZW)])) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + (k * cNXW) + (j * cNXW * cNZW)] +
						4 * cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						cVnp[i + ((k + 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == cNZV - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + ((k-1) * cNXW) + (j * cNXW * cNZW)] -
						4 * cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
						cVnp[i + ((k - 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + (k * cNXW) + (j * cNXW * cNZW)]
						+ (-4 * 0.025 + 1.0)*cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cVnp[i + ((k + 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(-2 * 0.025 / 3.0*cVnp[i + ((k + 3) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				else if (k == cNZV - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + ((k-1) * cNXW) + (j * cNXW * cNZW)]
						+ (4 * 0.09 + 1.0)*cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k-2) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cVnp[i + ((k - 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] +
						(2 * 0.09 / 3.0*cVnp[i + ((k - 3) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 4) * cNXW) + ((j)* cNXW * cNZW)]));
				}*/
				else
				{
					const double vijkp = cphiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];
					const double vijkpp = cphiInterpolateZ[i + ((k + 1)* cNXV) + (j* cNXV * (cNZV - 1))];

					const double vijkm = cphiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];
					const double vijkmm = cphiInterpolateZ[i + ((k - 2) * cNXV) + (j* cNXV * (cNZV - 1))];

					const double wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
					const double wijkpp = W[i + ((k+1)* cNXW) + (j*cNXW*cNZW)];

					const double wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];
					const double wijkmm = W[i + ((k - 2)* cNXW) + (j*cNXW*cNZW)];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*vijkmm * wijkmm 
						-Eimp / 2.0* vijkm * wijkm +
						Eimp / 2.0*vijkp * wijkp +
						Fimp / 4.0*vijkpp * wijkpp);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidxWConvSystemTD(
	const field3D& U,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXW * cNYW * cNZW;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYW; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = i + (k * cNXW) + (j * cNXW * cNZW);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == cNXW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == 1 || i == 2 || i == cNXW - 2 || i == cNXW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wipjk = cphiInterpolateX[i + (k * (cNXW-1)) + (j* (cNXW-1) * cNZW)];

					const double wimjk = cphiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const double uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const double uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					rhs[index] = ((wipjk * uipjk)
						- (wimjk * uimjk)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wipjk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i+1 + ((k)* cNXW) + ((j)* cNXW * cNZW)]
						);

					const double wimjk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i-1 + ((k) * cNXW) + ((j)* cNXW * cNZW)]
						);

					rhs[index] = ((wipjk * U[i + ((k)* cNXU) + ((j)* cNXU  * (cNZU-1))])
						- (wimjk * U[i-1 + ((k) * cNXU) + ((j)* cNXU * (cNZU - 1))])) /
						(ch);
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] +
						4 * cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)]);
				}
				else if (i == cNXW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						4 * cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] +
						cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZW)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)]
						+ (-4 * 0.025 + 1.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)] +
						(-2 * 0.025 / 3.0*cWnp[i + 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZW)]));
				}
				else if (i == cNXW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i-1 + (k * cNXU) + (j * cNXU * cNZW)]
						+ (4 * 0.09 + 1.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i -3 + (k * cNXU) + (j * cNXU * cNZW)] +
						(2 * 0.09 / 3.0*cWnp[i - 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZW)]));
				}*/
				else
				{
					const double wipjk = cphiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
					const double wippjk = cphiInterpolateX[i + 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const double wimjk = cphiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
					const double wimmjk = cphiInterpolateX[i - 2 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];

					const double uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const double uippjk = U[i+1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					const double uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
					const double uimmjk = U[i - 2 + (k*cNXU) + (j * cNXU * (cNZU - 1))];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*wimmjk * uimmjk -
						Eimp / 2.0*wimjk * uimjk +
						Eimp / 2.0*wipjk * uipjk +
						Fimp / 4.0*wippjk * uippjk);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidyWConvSystemTD(
	const field3D& V,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXW * cNYW * cNZW;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYW; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = i + (k * cNXW) + (j * cNXW * cNZW);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == cNYW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == 1 || j == 2 || j == cNYW - 2 || j == cNYW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wijpk = cphiInterpolateY[i + (j * cNXW) + (k * cNXW * cNZW)];

					const double wijmk = cphiInterpolateY[i + ((j - 1) * cNXW) + (k * cNXW * cNZW)];

					const double vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
					const double vijmk = V[i + (k * cNXV) + ((j-1) * cNXV * (cNZV - 1))];

					rhs[index] = ((wijpk * vijpk)
						- (wijmk * vijmk)) /
						(ch);
				}
				/*
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wijpk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k)* cNXW) + ((j+1)* cNXW * cNZW)]
						);

					const double wijmk = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k)* cNXW) + ((j-1)* cNXW * cNZW)]
						);

					rhs[index] = ((wijpk * V[i + ((k)* cNXV) + ((j)* cNXV * (cNZV - 1))])
						- (wijmk * V[i + ((k)* cNXV) + ((j-1)* cNXV * (cNZV-1))])) /
						(ch);
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + (j * cNXV * cNZW)] +
						4 * cWnp[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 1) * cNXV * cNZW)] -
						cWnp[i + (k * cNXW) + ((j + 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 2) * cNXV * cNZW)]);
				}
				else if (j == cNYW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]* V[i + (k * cNXV) + ((j-1) * cNXV * cNZW)] -
						4 * cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 2) * cNXV * cNZW)] +
						cWnp[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 3) * cNXV * cNZW)]);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 1) * cNXV * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + (j * cNXV * cNZW)]
						+ (-4 * 0.025 + 1.0)*cWnp[i + (k * cNXW) + ((j + 1)* cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 1)* cNXV * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cWnp[i + (k * cNXW) + ((j + 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 2) * cNXV * cNZW)] +
						(-2 * 0.025 / 3.0*cWnp[i + (k * cNXW) + ((j + 3) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 3) * cNXV * cNZW)]));
				}
				else if (j == cNYW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cWnp[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j) * cNXV * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + ((j-1) * cNXV * cNZW)]
						+ (4 * 0.09 + 1.0)*cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 2) * cNXV * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cWnp[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 3) * cNXV * cNZW)] +
						(2 * 0.09 / 3.0*cWnp[i + (k * cNXW) + ((j - 3) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 4) * cNXV * cNZW)]));
				}
				*/
				else
				{
					const double wijpk = cphiInterpolateY[i + (j * cNXW) + (k * cNXW * cNZW)];
					const double wijppk = cphiInterpolateY[i + ((j+1) * cNXW) + (k * cNXW * cNZW)];

					const double wijmk = cphiInterpolateY[i + ((j - 1) * cNXW) + (k * cNXW * cNZW)];
					const double wijmmk = cphiInterpolateY[i + ((j - 2) * cNXW) + (k * cNXW * cNZW)];

					const double vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
					const double vijppk = V[i + (k * cNXV) + ((j+1) * cNXV * (cNZV - 1))];

					const double vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];
					const double vijmmk = V[i + (k * cNXV) + ((j - 2) * cNXV * (cNZV - 1))];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*wijmmk * vijmmk -
						Eimp / 2.0*wijmk * vijmk +
						Eimp / 2.0*wijpk * vijpk +
						Fimp / 4.0*wijppk * vijppk);
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidzWConvSystemTD(
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = cNXW * cNYW * cNZW;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < cNYW; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = i + (k * cNXW) + (j * cNXW * cNZW);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == cNZW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == 1 || k == 2 || k == cNZW - 2 || k == cNZW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wijkp = cphiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];

					const double wijkm = cphiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];

					rhs[index] = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				/*else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wijkp = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k+1)* cNXW) + ((j)* cNXW * cNZW)]
						);

					const double wijkm = 0.5*(
						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k-1)* cNXW) + ((j)* cNXW * cNZW)]
						);

					rhs[index] = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] +
						4 * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == cNZW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] -
						4 * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] +
						cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 1.0 / ch * ((-2 * 0.025 / 3 - 1.0 / 3.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]
						+ (-4 * 0.025 + 1.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(-8 * 0.025 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(-2 * 0.025 / 3.0*cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				else if (k == cNZW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = -1.0 / ch * ((2 * 0.09 / 3 - 1.0 / 3.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNZW) + (j * cNXW * cNZW)]
						+ (4 * 0.09 + 1.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW *cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * 0.09 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(2 * 0.09 / 3.0*cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				*/
				else
				{
					const double wijkp = cphiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];
					const double wijkpp = cphiInterpolateZ[i + ((k+1) * cNXW) + (j * cNXW + (cNZW - 1))];

					const double wijkm = cphiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];
					const double wijkmm = cphiInterpolateZ[i + ((k - 2) * cNXW) + (j * cNXW + (cNZW - 1))];

					// im
					aa[index] = Dimp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = Dimp;

					rhs[index] = 1.0 / ch * (-Fimp / 4.0*wijkmm * wijkmm -
						Eimp / 2.0*wijkm * wijkm +
						Eimp / 2.0*wijkp * wijkp +
						Fimp / 4.0*wijkpp * wijkpp);
				}
			}

	return n2;
}

int FDColocated3DField::CreateInterpolateXSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (i == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX)+(j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + (k*NX)+(j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)])/2.0
								+bI * (phi[i - 1 + (k*NX) + (j*NX*NZ)] + phi[i + 2 + (k*NX) + (j*NX*NZ)])/2.0;
				}
			}

	return n2;
}

int FDColocated3DField::CreateInterpolateYSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j+1)*NX*NZ)]);
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j+1)*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + (k*NX) + ((j)*NX*NZ)] + phi[i + (k*NX) + ((j+1)*NX*NZ)]) / 2.0
						+ bI * (phi[i + (k*NX) + ((j-1)*NX*NZ)] + phi[i + (k*NX) + ((j+2)*NX*NZ)]) / 2.0;
				}
			}

	return n2;
}

int FDColocated3DField::CreateInterpolateZSystemTD(
	const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k+1)*NX) + (j*NX*NZ)]);
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k+1)*NX) + (j*NX*NZ)]);
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs[index] = aI * (phi[i + ((k)*NX) + (j*NX*NZ)] + phi[i + ((k+1)*NX) + (j*NX*NZ)]) / 2.0
						+ bI * (phi[i + ((k-1)*NX) + (j*NX*NZ)] + phi[i + ((k+2)*NX) + (j*NX*NZ)]) / 2.0;
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidxStaggeredSystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)])/ch;
				}
				else if (i == NXm - 1 || i == NXm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)])/ch;
				}
				else if (i == cNXP / 2 - 1  || i == cNXP/2 || i == cNXP / 2 + 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidyStaggeredSystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 || j == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j+1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)])/ch;
				}
				else if (j == NYm - 1 || j == NYm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j+1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)])/ch;
				}			
				else if (j == cNYP / 2 - 2 ||j == cNYP / 2 - 1 || j == cNYP / 2 || j == cNYP / 2 + 1 || j == cNYP / 2 + 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidzStaggeredSystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 || k == 1/* || true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)])/ch;
				}
				else if (k == NZm - 1 || k == NZm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)])/ch;
				}		
				else if (k == cNZP / 2 + 1 || k == cNZP/2 || k == cNZP / 2 - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + ((k+2)*NX) + (j*NX*NZ)] - phi[i + ((k - 1)*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k)*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidxStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const

{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = i + (k * NXm) + (j * NXm * NZ);

				if (i == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidyStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NYm = NY - 1;
	const int n2 = NX * NYm * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				if (j == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int FDColocated3DField::CreatedphidzStaggeredSystemTDSimplified(const field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZm);

				if (k == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs[index] = bII * (phi[i + ((k + 2)*NX) + (j*NX*NZ)] - phi[i + ((k - 1)*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k)*NX) + (j*NX*NZ)]) / ch;
				}
			}

	return n2;
}

int FDColocated3DField::Created2phidx2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i - 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + 1 + (k*NX) + (j*NX*NZ)])
						/(ch*ch);
				}
				else*/ if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be - phi at i+1;
						rhs[index] = (2.0 * -phi[i + 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i + 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					if (uvw != 0)
					{
						// force value to be - phi at i-1;
						rhs[index] = (2.0 * -phi[i - 1 + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i - 1 + (k*NX) + (j*NX*NZ)] +
							4.0 * phi[i - 2 + (k*NX) + (j*NX*NZ)] - phi[i - 3 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be - phi at 0;
						rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + 2 + (k*NX) + (j*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 2 + (k*NX) + (j*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i - 1 + (k*NX) + (j*NX*NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::Created2phidy2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + (k*NX) + ((j-1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + (k*NX) + ((j+1)*NX*NZ)])
						/ (ch*ch);
				}
				else */if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1)
					{
						rhs[index] = (2.0 * -phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j + 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j + 3)*NX*NZ)]) / (ch*ch);
					}
					
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1)
					{
						rhs[index] = (2.0 * -phi[i + (k*NX) + ((j - 1)*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + (k*NX) + ((j - 1)*NX*NZ)] +
							4.0 * phi[i + (k*NX) + ((j - 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 3)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					if (uvw != 1)
					{
						rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1)
					{
						rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (ch*ch);
					}
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + (k*NX) + ((j + 2)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 2)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + (k*NX) + ((j + 1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j - 1)*NX*NZ)]);
				}
			}

	return n2;
}

int FDColocated3DField::Created2phidz2SystemTD(const field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<double>& rhs)const
{
	const int n2 = NX * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = i + (k * NX) + (j * NX * NZ);

				/*if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs[index] = (phi[i + ((k-1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + ((k+1)*NX) + ((j)*NX*NZ)])
						/ (ch*ch);
				}
				else*/ if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 0)
					{
						// force mean value to be = ufreestream - LID IS MOVING
						const double phi0 = 2.0*cUfreestream - phi[i + ((k + 1)*NX) + ((j)*NX*NZ)];

						rhs[index] = (2.0 * phi0 - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else if (uvw == 1)
					{
						// zero mean value - lid moves only in x direction (u != 0);
						rhs[index] = (2.0 * -phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k + 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 2)
					{
						rhs[index] = (2.0 * -phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (2.0 * phi[i + (k*NX) + (j*NX*NZ)] - 5.0*phi[i + ((k - 1)*NX) + ((j)*NX*NZ)] +
							4.0 * phi[i + ((k - 2)*NX) + ((j)*NX*NZ)] - phi[i + ((k - 3)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 0)
					{
						// force mean value to be = ufreestream - LID IS MOVING
						const double phi0 = 2.0*cUfreestream - phi[i + ((k + 1)*NX) + ((j)*NX*NZ)];

						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							 + phi0) / (ch*ch);
					}
					else if (uvw == 1)
					{
						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							- phi[i + (k*NX) + (j*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 2)
					{
						rhs[index] = (-phi[i + (k*NX) + (j*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
					else
					{
						rhs[index] = (phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 *phi[i + (k*NX) + (j*NX*NZ)]
							+ phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]) / (ch*ch);
					}
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs[index] = (b / (4.0*ch*ch))*(phi[i + ((k + 2)*NX) + ((j)*NX*NZ)] - 2.0*phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 2)*NX) + ((j)*NX*NZ)])
						+ (a / (ch*ch))*(phi[i + ((k + 1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]);
				}
			}

	return n2;
}

void FDColocated3DField::EnforceContinuityVelnp(double dt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<double>& val,
	std::vector<double>& rhs, 
	LSSolver& lssolver)
{
	int nn = 0;
	if (true)
	{
		nn = CreatedphidxStaggeredSystemTDSimplified(cUnp, cNXU, cNYU, cNZU,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTDSimplified(cVnp, cNXV, cNYV, cNZV,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTDSimplified(cWnp, cNXW, cNYW, cNZW,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}
	else
	{
		nn = CreatedphidxStaggeredSystemTD(cUnp, cNXU, cNYU, cNZU,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTD(cVnp, cNXV, cNYV, cNZV,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTD(cWnp, cNXW, cNYW, cNZW,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}

	nn = CreateSparsePressureEquation(cdphidx, cdphidy, cdphidz,
		dt, ptr, col, val, rhs);

	if (cFirstIteration)
	{
		lssolver.PrecondtionCRS(nn, ptr, col, val);

		cFirstIteration = false;
	}

	double errOutPEqt = 0.0;
	int nitPEqt = 0;

	lssolver.SolvePreconditionedCRS(std::move(cPn), nn, ptr, col, val, rhs, errOutPEqt, nitPEqt);

	std::cout << "Pressure eq solved its " << nitPEqt << " err " << errOutPEqt << std::endl;
	   	

	// project p into uvw

	if (true)
	{
		nn = CreatedphidxStaggeredSystemTDSimplified(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTDSimplified(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTDSimplified(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}
	else
	{
		// set dp/dn = 0
		for (int j = 0; j < cNYP; j++)
			for (int k = 0; k < cNZP; k++)
			{
				cPn[0 + (k*cNXP) + (j*cNXP*cNZP)] = cPn[1 + (k*cNXP) + (j*cNXP*cNZP)];
				cPn[cNXP - 1 + (k*cNXP) + (j*cNXP*cNZP)] = cPn[cNXP - 2 + (k*cNXP) + (j*cNXP*cNZP)];
			}

		nn = CreatedphidxStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		// set dp/dn = 0
		for (int k = 0; k < cNZP; k++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn[i + (k*cNXP) + ((cNYP - 1)*cNXP*cNZP)] = cPn[i + (k*cNXP) + ((cNYP - 2)*cNXP*cNZP)];
				cPn[i + (k*cNXP) + (0*cNXP*cNZP)] = cPn[i + (k*cNXP) + (1*cNXP*cNZP)];
			}

		nn = CreatedphidyStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		// set dp/dn = 0
		for (int j = 0; j < cNYP; j++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn[i + (0*cNXP) + (j*cNXP*cNZP)] = cPn[i + (1*cNXP) + (j*cNXP*cNZP)];
				cPn[i + ((cNZP-1)*cNXP) + (j*cNXP*cNZP)] = cPn[i + ((cNZP - 2)*cNXP) + (j*cNXP*cNZP)];
			}

		nn = CreatedphidzStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}


	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				/*if (i == fdcolfield.ciCyInit - 1 || i == fdcolfield.ciCyEnd)
				{
					continue;
				}*/
				cUnp[i + (k*cNXU) + (j*cNXU*cNZU)] -=
					(dt / crho)*(cdphidx[i + (k*(cNXP - 1)) + (j*(cNXP - 1)*cNZP)]);
			}
		}
	}

	for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			for (int i = 1; i < cNXV - 1; i++)
			{
				/*if (j == fdcolfield.cjCyInit - 1 || j == fdcolfield.cjCyEnd)
				{
					continue;
				}*/
				cVnp[i + (k*cNXV) + (j*cNXV*cNZV)] -=
					(dt / crho)*(cdphidy[i + (k*cNXP) + ((j)*cNXP*cNZP)]);
			}
		}
	}

	for (int j = 1; j < cNYW - 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			for (int i = 1; i < cNXW - 1; i++)
			{
				/*if (k == fdcolfield.ckCyInit - 1 || k == fdcolfield.ckCyEnd)
				{
					continue;
				}*/
				cWnp[i + (k*cNXW) + (j*cNXW*cNZW)] -=
					(dt / crho)*(cdphidz[i + ((k)*cNXP) + (j*cNXP*(cNZP - 1))]);
			}
		}
	}
}