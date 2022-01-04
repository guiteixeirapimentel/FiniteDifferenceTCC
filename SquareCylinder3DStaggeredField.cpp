#include "SquareCylinder3DStaggeredField.h"
#include "LSSolver.h"

#include <string>
#include <thread>

SquareCylinder3DStaggeredField::SquareCylinder3DStaggeredField(double rho,
	double mu,
	double h,
	double l,
	double L,
	double H,
	double W,
	double ufreestream)
	:
	crho(rho),
	cmu(mu),
	ch(h),
	cL(L),
	cH(H),
	cW(W),

	cl(l),

	cFirstIteration(true),

	cDistRelFrenteRect(8),

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

	cNXCy(cl / ch),
	cNYCy(int(cNYP * 1.0)),
	cNZCy(cNXCy),

	ciCyInit(cDistRelFrenteRect* cNXCy),
	cjCyInit(/*(cNYP / 2 - cNYCy / 2)*/0),
	ckCyInit(cNZP / 2 - cNZCy / 2),

	ciCyEnd(cDistRelFrenteRect* cNXCy + cNXCy),
	cjCyEnd(/*(cNYP / 2 + cNYCy / 2) + (cNYCy % 2)*/cNYV),
	ckCyEnd((cNZP / 2 + cNZCy / 2) + (cNZCy % 2)),

	cUfreestream(ufreestream)
{
	cUn.resize(cNXU*cNYU*cNZU, cNXU,cNYU,cNZU, cUfreestream);
	cVn.resize(cNXV*cNYV*cNZV, cNXV,cNYV,cNZV, 0.0);
	cWn.resize(cNXW*cNYW*cNZW, cNXW,cNYW,cNZW, 0.0);
	cPn.resize(cNXP*cNYP*cNZP, cNXP,cNYP,cNZP, 0.0);

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
	//cPnp = cPn;

	cPhi = cPn;
	
	cdphidx = cPn;
	cdphidy.resize(cNXP*cNYP*cNZP, cNXP, cNYP, cNZP, 0.0);
	cdphidz.resize(cNXP*cNYP*cNZP, cNXP, cNYP, cNZP, 0.0);

	cd2phidx2.resize(cNXP*cNYP*cNZP, cNXP,cNYP,cNZP, 0.0);
	cd2phidy2.resize(cNXP*cNYP*cNZP, cNXP,cNYP,cNZP, 0.0);
	cd2phidz2.resize(cNXP*cNYP*cNZP, cNXP,cNYP,cNZP, 0.0);

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

SquareCylinder3DStaggeredField::~SquareCylinder3DStaggeredField()
{}

void SquareCylinder3DStaggeredField::SetBCFlatPlate()
{
	// freestream in all

	/*for (int j = 1; j < cNYU-1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			//cUn[0 + (k * cNXU) + (j*cNXU*cNZU)] = cUfreestream;
			cUn[1 + (k * cNXU) + (j*cNXU*cNZU)] = cUfreestream;
		}
	}*/

	for (int j = cjCyInit; j < cjCyEnd; j++)
		for (int k = ckCyInit; k < ckCyEnd; k++)
			for (int i = ciCyInit - 1; i < ciCyEnd; i++)
			{
				cUn(i, j, k) = 0.0;
			}

	for (int j = cjCyInit/* - 1*/; j < cjCyEnd; j++)
		for (int k = ckCyInit; k < ckCyEnd; k++)
			for (int i = ciCyInit; i < ciCyEnd; i++)
			{
				cVn(i, j, k) = 0.0;
			}

	for (int j = cjCyInit; j < cjCyEnd; j++)
		for (int k = ckCyInit - 1; k < ckCyEnd; k++)
			for (int i = ciCyInit; i < ciCyEnd; i++)
			{
				cWn(i, j, k) = 0.0;
			}

	cUnp = cUn;
	cVnp = cVn;
	cWnp = cWn;
}

void SquareCylinder3DStaggeredField::UpdateBCFlatPlate()
{
	// extrapolate outlet (i = NX - 1)

	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			{
				const int i = cNXU - 1;
				const int iu = cNXU - 2;

				cUnp(iu + 1, j, k) = cUnp(iu, j, k) +
					cVnp(i, j-1, k) - cVnp(i, j, k) +
					cWnp(i, j, k-1) - cWnp(i, j, k);
			}

			//cUnp[cNXU - 1 + (k * cNXU) + (j * cNXU * cNZU)] = cUnp[cNXU - 2 + (k * cNXU) + (j * cNXU * cNZU)];

		}
	}

	/*for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			cVn[cNXV - 1 + (k * cNXV) + (j * cNXV * cNZV)] = cVn[cNXV - 2 + (k * cNXV) + (j * cNXV * cNZV)];
		}
	}

	for (int j = 1; j < cNYW - 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			cWn[cNXW - 1 + (k * cNXW) + (j * cNXW * cNZW)] = cWn[cNXW - 2 + (k * cNXW) + (j * cNXW * cNZW)];
		}
	}*/


	/*for(int k = 1; k < cNZU - 1; k++)
		for (int i = 1; i < cNXU - 1; i++)
		{
			cUnp[i + (k * cNXU) + ((cNYU - 1) * cNXU * cNZU)] = cUnp[i + (k * cNXU) + ((cNYU - 2) * cNXU * cNZU)];

			cUnp[i + (k * cNXU) + (0 * cNXU * cNZU)] = cUnp[i + (k * cNXU) + (1 * cNXU * cNZU)];
		}


	for (int k = 1; k < cNZW - 1; k++)
		for (int i = 1; i < cNXW - 1; i++)
		{
			cWnp[i + (k * cNXW) + ((cNYW - 1) * cNXW * cNZW)] = cWnp[i + (k * cNXW) + ((cNYW - 2) * cNXW * cNZW)];
			cWnp[i + (k * cNXW) + (0 * cNXW * cNZW)] = cWnp[i + (k * cNXW) + (1 * cNXW * cNZW)];
		}



	for (int j = 1; j < cNYU - 1; j++)
		for (int i = 1; i < cNXU - 1; i++)
		{
			cUnp[i + ((cNZU - 1) * cNXU) + (j * cNXU * cNZU)] = cUnp[i + ((cNZU - 2) * cNXU) + (j * cNXU * cNZU)];
			cUnp[i + (0 * cNXU) + (j * cNXU * cNZU)] = cUnp[i + (1 * cNXU) + (j * cNXU * cNZU)];
		}

	for (int j = 1; j < cNYV - 1; j++)
		for (int i = 1; i < cNXV - 1; i++)
		{
			cVnp[i + ((cNZV - 1) * cNXW) + (j * cNXW * cNZW)] = cVnp[i + ((cNZV - 2) * cNXW) + (j * cNXW * cNZW)];
			cVnp[i + (0 * cNXW) + (j * cNXW * cNZW)] = cVnp[i + (1 * cNXW) + (j * cNXW * cNZW)];
		}*/
}

int SquareCylinder3DStaggeredField::CreatedphidxSystemTDOUCS3(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
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
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0 && false)
						rhs(i, j, k) = 0.5 / ch * (-3 * -phi(i + 1, j, k) +
							4 * phi(i + 1, j, k) - phi(i + 2, j, k));
					else
						rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
							4 * phi(i + 1, j, k) - phi(i + 2, j, k));
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0 && false)
						rhs(i, j, k) = 0.5 / ch * (3 * -phi(i - 1, j, k) -
							4 * phi(i - 1, j, k) + phi(i - 2, j, k));
					else
						rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
							4 * phi(i - 1, j, k) + phi(i - 2, j, k));
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0 && false)
						rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betaj2 + 1.0)*phi(i + 1, j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2, j, k) + (2 * betaj2 / 3.0*phi(i + 3, j, k)));
					else
						rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i - 1, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betaj2 + 1.0)*phi(i + 1, j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2, j, k) + (2 * betaj2 / 3.0*phi(i + 3,j, k)));
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0 && false)
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k)+ (2 * betajm2 / 3.0*phi(i - 3, j, k)));
					else
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k)+ (2 * betajm2 / 3.0*phi(i - 3, j, k)));
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k) + (2 * betajm2 / 3.0*phi(i - 3, j, k)));
				}
				else if (uvw == 0 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i - 1, j, k) + phi(i - 2, j, k));
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i - 1, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i + 1, j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2, j, k) + (2 * betaj2 / 3.0*phi(i + 3, j, k)));
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i + 1, j, k) - phi(i + 2, j, k));
				}
				else if (uvw == 1 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k) + (2 * betajm2 / 3.0*phi(i - 3, j, k)));
				}
				else if (uvw == 1 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i - 1, j, k) + phi(i - 2, j, k));
				}
				else if (uvw == 1 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i - 1 , j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i + 1 , j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2 , j, k) + (2 * betaj2 / 3.0*phi(i + 3 , j, k)));
				}
				else if (uvw == 1 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i + 1 , j, k) - phi(i + 2 , j, k));
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1 , j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i - 1 , j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2 , j, k) + (2 * betajm2 / 3.0*phi(i - 3 , j, k)));
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i - 1 , j, k) + phi(i - 2 , j, k));
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i - 1 , j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i + 1 , j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2 , j, k) + (2 * betaj2 / 3.0*phi(i + 3 , j, k)));
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i + 1 , j, k) - phi(i + 2 , j, k));
				}
				else
				{
					// im
					aa[index] = pim;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = pip;

					rhs(i, j, k) = 1.0 / ch *
						(qm2*phi(i - 2 , j, k) + qm1 * phi(i - 1 , j, k) +
							+q0 * phi(i, j, k) +
							qp1 * phi(i + 1 , j, k) + qp2 * phi(i + 2 , j, k));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidySystemTDOUCS3(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs)const
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
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1 && false)
						rhs(i, j, k) = 0.5 / ch * (-3 * -phi(i, ((j + 1)), k) +
							4 * phi(i, j + 1, k) - phi(i, j + 2, k));
					else
						rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
							4 * phi(i, j + 1, k) - phi(i, j+2, k));
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1 && false)
						rhs(i, j, k) = 0.5 / ch * (3 * -phi(i, j - 1, k) -
							4 * phi(i, j - 1, k) + phi(i, j - 2, k));
					else
						rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
							4 * phi(i, j - 1, k) + phi(i, j - 2, k));
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1 && false)
						rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betaj2 + 1.0)*phi(i, j + 1, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, ((j + 2))) + (2 * betaj2 / 3.0*phi(i, j, j + 3)));
					else
						rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j - 1, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betaj2 + 1.0)*phi(i, j + 1, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j + 2, k) + (2 * betaj2 / 3.0*phi(i, j + 3, k)));
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 1 && false)
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i, j - 1, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j - 2, k) + (2 * betajm2 / 3.0*phi(i, j - 3, k)));
					else
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j + 1, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i, j - 1, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j - 2, k) + (2 * betajm2 / 3.0*phi(i, j - 3, k)));
				}
				else if (uvw == 0 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j + 1, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j - 1, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j - 2, k) + (2 * betajm2 / 3.0*phi(i, j - 3, k)));
				}
				else if (uvw == 0 && j == cjCyInit && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i, j - 1, k) + phi(i, j - 2, k));
				}
				else if (uvw == 0 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j - 1, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j + 1, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j + 2, k) + (2 * betaj2 / 3.0*phi(i, j + 3, k)));
				}
				else if (uvw == 0 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i, j + 1, k) - phi(i, j + 2, k));
				}
				else if (uvw == 1 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j + 1, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j - 1, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j - 2, k) + (2 * betajm2 / 3.0*phi(i, j - 3, k)));
				}
				else if (uvw == 1 && j == cjCyInit && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i, j - 1, k) + phi(i, j - 2, k));
				}
				else if (uvw == 1 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j - 1, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j + 1, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j + 2, k) + (2 * betaj2 / 3.0*phi(i, j + 3, k)));
				}
				else if (uvw == 1 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i, j + 1, k) - phi(i, j + 2, k));
				}
				else if (uvw == 2 && j == cjCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j + 1, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j - 1, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j - 2, k) + (2 * betajm2 / 3.0*phi(i, j - 3, k)));
				}
				else if (uvw == 2 && j == cjCyInit && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i, j - 1, k) + phi(i, j - 2, k));
				}
				else if (uvw == 2 && j == cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j - 1, k)- (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j + 1, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j + 2, k) + (2 * betaj2 / 3.0*phi(i, j + 3, k)));
				}
				else if (uvw == 2 && j == cjCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i, j + 1, k) - phi(i, j + 2, k));
				}
				else
				{
					// im
					aa[index] = pim;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = pip;

					rhs(i, j, k) = 1.0 / ch *
						(qm2*phi(i, j - 2, k) + qm1 * phi(i, j - 1, k) +
							q0 * phi(i, j, k) +
							qp1 * phi(i, j + 1, k) + qp2 * phi(i, j + 2, k));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidzSystemTDOUCS3(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
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
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2 && false)
						rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
							4 * phi(i, j, ((k + 1))) - phi(i, j, (k + 2)));
					else
						rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
							4 * phi(i, j, k + 1) - phi(i, j, (k + 2)));
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2 && false)
						rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) - 4 * phi(i, j, k - 1) + phi(i, j, k - 2));
					else
						rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) - 4 * phi(i, j, (k - 1)) + phi(i, j, k - 2));
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 2 && false)
						rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betaj2 + 1.0)*phi(i, j, k + 1) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, k + 2) + (2 * betaj2 / 3.0*phi(i, j, k + 3)));
					else
						rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j, k - 1) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betaj2 + 1.0)*phi(i, j, k + 1) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, k + 2) + (2 * betaj2 / 3.0*phi(i, j, k + 3)));
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;
					if (uvw != 2 && false)
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
					else
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j, k + 1) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
				}
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j, k + 1) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i, j, k - 1) + phi(i, j, k - 2));
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i, j, k + 1) - phi(i,j, k + 2));
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j, k - 1) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j, k + 1) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, k + 2) + (2 * betaj2 / 3.0*phi(i, j, k + 3)));
				}
				else if (uvw == 1 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j, k + 1) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
				}
				else if (uvw == 1 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i, j, k - 1) + phi(i, j, k - 2));
				}
				else if (uvw == 1 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i, j, k + 1) - phi(i, j, k + 2));
				}
				else if (uvw == 1 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j, k - 1) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j, k + 1) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, k + 2) + (2 * betaj2 / 3.0*phi(i, j, k + 3)));
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j, k + 1) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
				}
				else if (uvw == 2 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) - 4 * phi(i, j, k - 1) + phi(i, j, k - 2));
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) + 4 * phi(i, j, k + 1) - phi(i, j, k + 2));
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j, k - 1) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j, k + 1) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, k + 2) + (2 * betaj2 / 3.0*phi(i, j, k + 3)));
				}
				else
				{
					// im
					aa[index] = pim;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = pip;

					rhs(i, j, k) = 1.0 / ch *
						(qm2*phi(i, j, k - 2) + qm1 * phi(i, j, k - 1)
							+ q0 * phi(i, j, k)
							+ qp1 * phi(i, j, k + 1) + qp2 * phi(i, j, k + 2));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidxSystemTDSUCS(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NYc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				double velConv = 0.0;

				if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1 || k == 0 || k == NZ - 1)
				{
					velConv = NAN;
				}
				else if (uvw == 0)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j - 1, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k - 1);
					}
				}
				else if (uvw == 1)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i - 1, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k - 1);
					}
				}
				else if (uvw == 2)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i - 1, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j - 1, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k);
					}
				}

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));*/

					rhs(i, j, k) = (phi(i+1, j, k) - phi(i, j, k))/ch;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));*/

					rhs(i, j, k) = (phi(i, j, k) - phi(i-1, j, k)) / ch;
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i-1, j, k)) / (ch/2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / (ch / 2.0);
					}
				}
				else if (uvw == 0 && i >= ciCyInit && i < ciCyEnd - 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 1 && i >= ciCyInit + 1  && i < ciCyEnd - 1 && j >= cjCyInit - 1 && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit + 1 && i < ciCyEnd - 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));*/

					rhs(i, j, k) = (0.0 - phi(i - 1, j, k)) / ch;
				}
				else if (uvw == 0 && i == ciCyInit - 2 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));*/

					rhs(i, j, k) = (phi(i + 1, j, k) - 0.0) / ch;
				}
				else if (uvw == 1 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i -1, j, k)
						- 6 * phi(i -2, j, k) + phi(i-3, j, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					}
					else
					{
						// derivada calculada entre parede e freestream
						rhs(i, j, k) = (/*phi(i+1, j, k)*/0.0 - phi(i, j, k)) / (ch/2.0);
					}
				}
				else if (uvw == 1 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i-1, j, k)
						+ 48 * phi(i -1, j, k) - 36 * phi(i-2, j, k)
						+ 16 * phi(i-3,j,k) - 3 * phi(i-4, j, k));*/
					rhs(i, j, k) = (/*phi(i, j, k)*/0.0 - phi(i-1, j, k)) / (ch/2);
				}
				else if (uvw == 1 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i+1,j,k)
						- 6 * phi(i+2, j, k) + phi(i + 3,j,k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0/*phi(i - 1, j, k)*/) / (ch/2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 1 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i + 1, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));*/

					rhs(i, j, k) = (phi(i+1, j, k) - 0.0/*phi(i, j, k)*/) / (ch/2.0);
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));*/
					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (/*phi(i + 1, j, k)*/0.0 - phi(i, j, k)) / (ch/2);
					}
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i-1, j, k)
						+ 48 * phi(i-1, j, k) - 36 * phi(i-2, j, k)
						+ 16 * phi(i-3, j, k) - 3 * phi(i-4, j, k));*/

					rhs(i, j, k) = (0.0 - phi(i - 1, j, k)) / (ch/2.0);
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i+1, j, k)
						- 6 * phi(i+2, j, k) + phi(i+3, j, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch/2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i+1, j, k)
						+ 48 * phi(i+1, j, k) - 36 * phi(i+2, j, k)
						+ 16 * phi(i+3, j, k) - 3 * phi(i+4, j, k));*/

					rhs(i, j, k) = (phi(i+1, j, k) - 0.0) / (ch / 2.0);
				}
				else
				{
					if (velConv >= 0.0)
					{
						// im
						aa[index] = velPosSUCSalfam1;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = velPosSUCSalfap1;

						rhs(i, j, k) = 1.0 / ch *
							(velPosSUCScm2 * phi(i-2, j, k)
								+ velPosSUCScm1 * phi(i-1, j, k)
								+ velPosSUCSc0 * phi(i, j, k)
								+ velPosSUCSc1 * phi(i+1, j, k));

						//// im
						//aa[index] = 0.0;

						//// i
						//bb[index] = 1.0;

						//// ip
						//cc[index] = 0.0;

						//rhs(i, j, k) =  (phi[i + (k*NX) + (j*NX*NZ)]-
						//		phi[i - 1 + (k*NX) + (j*NX*NZ)])/ch
						//		;
					}
					else
					{
						// im
						aa[index] = velNegSUCSalfam1;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = velNegSUCSalfap1;

						rhs(i, j, k) = 1.0 / ch *
							(velNegSUCSc2 * phi(i+2, j, k)
								+ velNegSUCSc1 * phi(i+1, j, k)
								+ velNegSUCSc0 * phi(i, j, k)
								+ velNegSUCScm1 * phi(i-1, j, k)
								);

						//// im
						//aa[index] = 0.0;

						//// i
						//bb[index] = 1.0;

						//// ip
						//cc[index] = 0.0;

						//rhs(i, j, k) = (phi[i + 1 + (k*NX) + (j*NX*NZ)] -
						//	phi[i + (k*NX) + (j*NX*NZ)]) / ch
						//	;
					}

				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidySystemTDSUCS(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NYc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				double velConv;

				if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1 || k == 0 || k == NZ - 1)
				{
					velConv = NAN;
				}
				else if (uvw == 0)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j - 1, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k - 1);
					}
				}
				else if (uvw == 1)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i - 1, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k - 1);
					}
				}
				else if (uvw == 2)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i - 1, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j - 1, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k);
					}
				}

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j+1, k) - 36 * phi(i, j+2, k)
						+ 16 * phi(i, j+3, k) - 3 * phi(i, j+4, k));*/

					rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / ch;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j-2, k)
						+ 16 * phi(i, j-3, k) - 3 * phi(i, j-4, k));*/

					rhs(i, j, k) = (phi(i, j, k) - phi(i, j-1, k)) / ch;
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j-1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j+1, k)
						- 6 * phi(i, j+2, k) + phi(i, j+3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j+1, k) - phi(i, j, k)) / ch;
					}
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j+1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j-1, k)
						- 6 * phi(i, j-2, k) + phi(i, j-3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j - 1, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}
				}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit + 1 && j < cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit + 1 && j < cjCyEnd -1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit + 1 && j < cjCyEnd - 1 && k >= ckCyInit-1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j+1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j-1, k)
						- 6 * phi(i, j-2, k) + phi(i, j-3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j - 1, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch/2.0);
					}
				}
				else if (uvw == 0 && j == cjCyInit && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j-1, k) - 36 * phi(i, j-2, k)
						+ 16 * phi(i, j-3, k) - 3 * phi(i, j-4, k));*/

					rhs(i, j, k) = (0.0 - phi(i, j - 1, k)) / (ch / 2.0);
				}
				else if (uvw == 0 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j-1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j+1, k)
						- 6 * phi(i, j+2, k) + phi(i, j+3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 0 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j-1, k)
						+ 48 * phi(i, j+1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));*/

					rhs(i, j, k) = (phi(i+1, j, k) - 0.0) / (ch / 2.0);
				}
				else if (uvw == 1 && j == cjCyInit - 2 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j + 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j - 1, k)
						- 6 * phi(i, j - 2, k) + phi(i, j - 3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j - 1, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 1 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j - 2, k)
						+ 16 * phi(i, j - 3, k) - 3 * phi(i, j -4, k));*/

					rhs(i, j, k) = (0.0 - phi(i, j-1, k)) / ch;
				}
				else if (uvw == 1 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j - 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j + 1, k)
						- 6 * phi(i, j + 2, k) + phi(i, j + 3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 1 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j + 1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));*/

					rhs(i, j, k) = (phi(i, j + 1, k) - 0.0) / ch;
				}
				else if (uvw == 2 && j == cjCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j - 1, k)
						- 6 * phi(i, j - 2, k) + phi(i, j -3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j - 1, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}
				}
				else if (uvw == 2 && j == cjCyInit && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i, j - 1, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j - 2, k)
						+ 16 * phi(i, j - 3, k) - 3 * phi(i, j - 4, k));*/

					rhs(i, j, k) = (0.0 - phi(i, j - 1, k)) / (ch / 2.0);

				}
				else if (uvw == 2 && j == cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j + 1, k)
						- 6 * phi(i, j +2, k) + phi(i, j + 3, k));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j+1, k) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 2 && j == cjCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i, j + 1, k)
						+ 48 * phi(i, j + 1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));*/

					rhs(i, j, k) = (phi(i, j+1, k) - 0.0) / (ch / 2.0);
				}
				else
				{
					if (velConv > 0.0)
					{
						// im
						aa[index] = velPosSUCSalfam1;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = velPosSUCSalfap1;

						rhs(i, j, k) = 1.0 / ch *
							(velPosSUCScm2 * phi(i, j-2, k)
								+ velPosSUCScm1 * phi(i, j-1, k)
								+ velPosSUCSc0 * phi(i, j, k)
								+ velPosSUCSc1 * phi(i, j+1, k));

						//// im
						//aa[index] = 0.0;

						//// i
						//bb[index] = 1.0;

						//// ip
						//cc[index] = 0.0;

						//rhs(i, j, k) = 1.0 / ch *
						//	(-phi[i + (k*NX) + ((j - 1)*NX*NZ)]
						//		+  phi[i + (k*NX) + (j*NX*NZ)]);
					}
					else
					{
						// im
						aa[index] = velNegSUCSalfam1;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = velNegSUCSalfap1;

						rhs(i, j, k) = 1.0 / ch *
							(velNegSUCSc2 * phi(i, j+2, k)
								+ velNegSUCSc1 * phi(i, j+1, k)
								+ velNegSUCSc0 * phi(i, j, k)
								+ velNegSUCScm1 * phi(i, j-1, k)
								);

						//// im
						//aa[index] = 0.0;

						//// i
						//bb[index] = 1.0;

						//// ip
						//cc[index] = 0.0;

						//rhs(i, j, k) = 1.0 / ch *
						//	(phi[i + (k*NX) + ((j + 1)*NX*NZ)]
						//		- phi[i + (k*NX) + (j*NX*NZ)]);
					}
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidzSystemTDSUCS(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NYc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				double velConv;

				if (i == 0 || i == NX - 1 || j == 0 || j == NY - 1 || k == 0 || k == NZ - 1)
				{
					velConv = NAN;
				}
				else if (uvw == 0)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j - 1, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k - 1);
					}
				}
				else if (uvw == 1)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i - 1, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k - 1);
					}
				}
				else if (uvw == 2)
				{
					if (uvwc == 0)
					{
						velConv = convVel(i - 1, j, k);
					}
					else if (uvwc == 1)
					{
						velConv = convVel(i, j - 1, k);
					}
					else if (uvwc == 2)
					{
						velConv = convVel(i, j, k);
					}
				}

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k+1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));*/

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));*/

					rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k-1)) / ch;
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j, k - 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k+1) - phi(i, j, k)) / ch;
					}
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}
				}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));*/

					rhs(i, j, k) = (0.0 - phi(i, j, k-1)) / (ch / 2.0);
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));*/

					rhs(i, j, k) = (phi(i, j, k+1) - 0.0) / (ch / 2.0);
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 1 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k -2) + phi(i, j, k - 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}
				}
				else if (uvw == 1 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));*/

					rhs(i, j, k) = (0.0 - phi(i, j, k - 1)) / (ch / 2.0);
				}
				else if (uvw == 1 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));*/

					rhs(i, j, k) = (phi(i, j, k + 1) - 0.0) / (ch / 2.0);
				}
				else if (uvw == 1 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 2 && k == ckCyInit - 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / ch;
					}
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));*/

					rhs(i, j, k) = (0.0 - phi(i, j, k - 1)) / ch;
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));*/

					rhs(i, j, k) = (phi(i, j, k+1) - 0.0) / ch;
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					/*rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j, k - 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));*/

					if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
					}
				}
				else
				{
				
					if (velConv >= 0.0)
					{
						// im
						aa[index] = velPosSUCSalfam1;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = velPosSUCSalfap1;

						rhs(i, j, k) = 1.0 / ch *
							(velPosSUCScm2 * phi(i, j, k-2)
								+ velPosSUCScm1 * phi(i, j, k-1)
								+ velPosSUCSc0 * phi(i, j, k)
								+ velPosSUCSc1 * phi(i, j, k+1));

						//// im
						//aa[index] = 0.0;

						//// i
						//bb[index] = 1.0;

						//// ip
						//cc[index] = 0.0;

						//rhs(i, j, k) = 1.0 / ch *
						//	(-phi[i + ((k - 1)*NX) + ((j)*NX*NZ)]
						//		+ phi[i + (k*NX) + (j*NX*NZ)]);
					}
					else
					{
						// im
						aa[index] = velNegSUCSalfam1;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = velNegSUCSalfap1;

						rhs(i, j, k) = 1.0 / ch *
							(velNegSUCSc2 * phi(i, j, k+2)
								+ velNegSUCSc1 * phi(i, j, k+1)
								+ velNegSUCSc0 * phi(i, j, k)	
								+ velNegSUCScm1 * phi(i, j, k-1)
								);

						//// im
						//aa[index] = 0.0;

						//// i
						//bb[index] = 1.0;

						//// ip
						//cc[index] = 0.0;

						//rhs(i, j, k) = 1.0 / ch *
						//	(phi[i + ((k + 1)*NX) + ((j)*NX*NZ)]
						//		- phi[i + (k*NX) + (j*NX*NZ)]);
					}
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidxSystemTDUPWIND(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NYc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);
				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 0 && i >= ciCyInit && i < ciCyEnd - 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 1 && i >= ciCyInit + 1 && i < ciCyEnd - 1 && j >= cjCyInit - 1 && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit + 1 && i < ciCyEnd - 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (uvw == 0 && i == ciCyInit - 2 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));
				}
				else if (uvw == 1 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 1 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i - 1, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (uvw == 1 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (uvw == 1 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i + 1, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i - 1, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i + 1, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));
				}
				else
				{
					double velConv;

					if (uvw == 0)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j - 1, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k - 1);
						}
					}
					else if (uvw == 1)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i - 1, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k - 1);
						}
					}
					else if (uvw == 2)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i - 1, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j - 1, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k);
						}
					}

					if (velConv >= 0.0)
					{
						// im
						aa[index] = 0.0;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = 0.0;

						rhs(i, j, k) = 1.0 / ch *
							(
								- phi(i - 1, j, k)
								+ phi(i, j, k)
							);
					}
					else
					{
						// im
						aa[index] = 0.0;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = 0.0;

						rhs(i, j, k) = 1.0 / ch *
							(
								phi(i + 1, j, k)
								- phi(i, j, k)
								);

					}

				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidySystemTDUPWIND(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NYc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j + 1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j - 2, k)
						+ 16 * phi(i, j - 3, k) - 3 * phi(i, j - 4, k));
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i, j - 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j + 1, k)
						- 6 * phi(i, j + 2, k) + phi(i, j + 3, k));
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i, j + 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j - 1, k)
						- 6 * phi(i, j - 2, k) + phi(i, j - 3, k));
				}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit + 1 && j < cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit + 1 && j < cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit + 1 && j < cjCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i, j + 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j - 1, k)
						- 6 * phi(i, j - 2, k) + phi(i, j - 3, k));
				}
				else if (uvw == 0 && j == cjCyInit && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j - 2, k)
						+ 16 * phi(i, j - 3, k) - 3 * phi(i, j - 4, k));
				}
				else if (uvw == 0 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i, j - 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j + 1, k)
						- 6 * phi(i, j + 2, k) + phi(i, j + 3, k));
				}
				else if (uvw == 0 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j - 1, k)
						+ 48 * phi(i, j + 1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));
				}
				else if (uvw == 1 && j == cjCyInit - 2 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i, j + 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j - 1, k)
						- 6 * phi(i, j - 2, k) + phi(i, j - 3, k));
				}
				else if (uvw == 1 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j - 2, k)
						+ 16 * phi(i, j - 3, k) - 3 * phi(i, j - 4, k));
				}
				else if (uvw == 1 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i, j - 1, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j + 1, k)
						- 6 * phi(i, j + 2, k) + phi(i, j + 3, k));
				}
				else if (uvw == 1 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j + 1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));
				}
				else if (uvw == 2 && j == cjCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j - 1, k)
						- 6 * phi(i, j - 2, k) + phi(i, j - 3, k));
				}
				else if (uvw == 2 && j == cjCyInit && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i, j - 1, k)
						+ 48 * phi(i, j - 1, k) - 36 * phi(i, j - 2, k)
						+ 16 * phi(i, j - 3, k) - 3 * phi(i, j - 4, k));
				}
				else if (uvw == 2 && j == cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j + 1, k)
						- 6 * phi(i, j + 2, k) + phi(i, j + 3, k));
				}
				else if (uvw == 2 && j == cjCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i, j + 1, k)
						+ 48 * phi(i, j + 1, k) - 36 * phi(i, j + 2, k)
						+ 16 * phi(i, j + 3, k) - 3 * phi(i, j + 4, k));
				}
				else
				{
					double velConv;

					if (uvw == 0)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j - 1, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k - 1);
						}
					}
					else if (uvw == 1)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i - 1, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k - 1);
						}
					}
					else if (uvw == 2)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i - 1, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j - 1, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k);
						}
					}

					if (velConv > 0.0)
					{
						// im
						aa[index] = 0.0;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = 0.0;

						rhs(i, j, k) = 1.0 / ch *
							(
								-phi(i, j - 1, k)
								+ phi(i, j, k)
							);

					}
					else
					{
						// im
						aa[index] = 0.0;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = 0.0;

						rhs(i, j, k) = 1.0 / ch *
							(
								phi(i, j + 1, k)
								- phi(i, j, k)
							);					
					}
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidzSystemTDUPWIND(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NYc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i, j, k - 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));
				}
				else if (uvw == 1 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 1 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (uvw == 1 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (uvw == 1 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));
				}
				else if (uvw == 2 && k == ckCyInit - 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i, j, k - 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));
				}
				else
				{
					double velConv;

					if (uvw == 0)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j - 1, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k - 1);
						}
					}
					else if (uvw == 1)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i - 1, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k - 1);
						}
					}
					else if (uvw == 2)
					{
						if (uvwc == 0)
						{
							velConv = convVel(i - 1, j, k);
						}
						else if (uvwc == 1)
						{
							velConv = convVel(i, j - 1, k);
						}
						else if (uvwc == 2)
						{
							velConv = convVel(i, j, k);
						}
					}

					if (velConv >= 0.0)
					{
						// im
						aa[index] = 0.0;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = 0.0;

						rhs(i, j, k) = 1.0 / ch *
							(
								-phi(i, j, k - 1)
								+phi(i, j, k)
							);
						
					}
					else
					{
						// im
						aa[index] = 0.0;

						// i
						bb[index] = 1.0;

						// ip
						cc[index] = 0.0;

						rhs(i, j, k) = 1.0 / ch *
							(
								+phi(i, j, k + 1)
								-phi(i, j, k)
								);						
					}
				}
			}

	return n2;
}


int SquareCylinder3DStaggeredField::Created2phidx2SystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				/*if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi[i - 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + 1 + (k*NX) + (j*NX*NZ)])
						/(ch*ch);
				}
				else*/
				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i+1, j, k) +
						4.0 * phi(i + 2, j, k) - phi(i + 3, j, k)) / (ch*ch);
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i - 1, j, k) +
						4.0 * phi(i - 2, j, k) - phi(i - 3, j, k)) / (ch*ch);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
						+ phi(i - 1, j, k)) / (ch*ch);
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
						+ phi(i - 1, j, k)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i + 2, j, k) - 2.0*phi(i, j, k) + phi(i -2, j, k))
						+ (a / (ch*ch))*(phi(i + 1, j, k) - 2.0 * phi(i, j, k) + phi(i - 1, j, k));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::Created2phidy2SystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				/*if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi[i + (k*NX) + ((j-1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + (k*NX) + ((j+1)*NX*NZ)])
						/ (ch*ch);
				}
				else */if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j + 1, k) +
						4.0 * phi(i, j + 2, k) - phi(i, j + 3, k)) / (ch*ch);
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j - 1, k) +
						4.0 * phi(i, j - 2, k) - phi(i, j - 3, k)) / (ch*ch);
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - 2.0 *phi(i, j, k)
						+ phi(i, j - 1, k)) / (ch*ch);
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - 2.0 *phi(i, j, k)
						+ phi(i, j - 1, k)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i, j + 2, k) - 2.0*phi(i, j, k) + phi(i, j - 2, k))
						+ (a / (ch*ch))*(phi(i, j + 1, k) - 2.0 * phi(i, j, k) + phi(i, j - 1, k));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::Created2phidz2SystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * NY * NZ;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				/*if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi[i + ((k-1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + ((k+1)*NX) + ((j)*NX*NZ)])
						/ (ch*ch);
				}
				else*/ if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j, k + 1) +
						4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j, k - 1) +
						4.0 * phi(i, j, k -2) - phi(i, j, k - 3)) / (ch*ch);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i, j, k + 2) - 2.0*phi(i, j, k) + phi(i, j, k -2))
						+ (a / (ch*ch))*(phi(i, j, k + 1) - 2.0 * phi(i, j, k) + phi(i, j, k - 1));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreateSparsePressureEquation(const Field3D& dUdx,
	const Field3D& dVdy,
	const Field3D& dWdz,
	const double dt,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<double>& val,
	Field3DIKJ& rhs)const
{
	const int n2 = cNXP * cNYP * cNZP; // Number of points in the grid.

	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 7);

	val.clear();
	val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2, cNXP, cNYP, cNZP, 0.0);

	for (int j = 0; j < cNYP; ++j)
	{
		for (int k = 0; k < cNZP; k++)
		{
			for (int i = 0; i < cNXP; i++)
			{
				if (i == cNXP - 2)
				{
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
				}
				//else if (
				//	// baixo
				//	(i == 0 && k == 0) || (i == cNXP - 1 && k == 0)
				//	|| (j == cNYP - 1 && k == 0) || (k == 0 && j == 0) ||
				//	// cima
				//	(i == 0 && k == cNZP - 1) || (i == cNXP - 1 && k == cNZP - 1)
				//	|| (j == cNYP - 1 && k == cNZP - 1) || (k == cNZP - 1 && j == 0) ||
				//	//laterais
				//	(i == 0 && j == 0) || (i == cNXP - 1 && j == 0)
				//	|| (i == 0 && j == cNYP - 1) || (i == cNXP - 1 && j == cNYP - 1)
				//	)
				//{
				//	col.push_back(index);
				//	val.push_back(1.0);
				//	rhs(i, j, k) = 0.0;
				//}
				else if (i == cNXP - 1)
				{
					// Outlet => Pressure set to zero
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
				}
				else if (k == 0)
				{
					// surface => dPdn = 0
					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
				}
				else if (k == cNZP - 1)
				{
					// top => dPdn = 0

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
				}
				else if (i == 0)
				{
					// inlet => dPdn = 0

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;

				}
				else if (j == cNYP - 1)
				{
					// back => dPdn = 0

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
				}
				else if (j == 0)
				{
					// front => dPdn = 0

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
				}
				else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					// dentro do sólido

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(1.0);

					rhs(i, j, k) = 0.0;
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
					col.push_back(rhs.GetIndex(i, j, k - 1));
					if (k == 1 || (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkpi
					col.push_back(rhs.GetIndex(i, j, k+1));
					if (k == cNZP - 2 || (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkim
					col.push_back(rhs.GetIndex(i-1, j, k));
					if (i == 1 || (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}


					// jkip
					col.push_back(rhs.GetIndex(i+1, j, k));
					if (i == cNXP - 2 || (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

					// jmki
					col.push_back(rhs.GetIndex(i, j-1, k));
					if (j == 1 || (j == cjCyEnd && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}


					// jpki
					col.push_back(rhs.GetIndex(i, j+1, k));
					if (j == cNYP - 2 || (j == cjCyInit - 1 && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(valijk);

					rhs(i, j, k) = ch * ch*(crho / dt)*(
						dUdx(i-1, j, k) +
						dVdy(i, j-1, k) +
						dWdz(i, j, k-1)
						);
				}
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

					rhs(i, j, k) = ch * ch*(crho / dt)*(
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

//void SquareCylinder3DStaggeredField::CalcdphidtUmomDiv(field3D& dphidt,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc)
//{
//	int nn = 0;
//
//	static field3D vnpInterpolated;
//	static field3D wnpInterpolated;
//	static field3D unpXinterpolated;
//	static Field3DJIK unpYinterpolated;
//	static Field3DKJI unpZinterpolated;
//
//	static field3D dphidx, dphidy, dphidz;
//	static Field3DIKJ d2phidx2;
//	static Field3DJIK d2phidy2; 
//	static Field3DKJI d2phidz2;
//
//	/*
//		vnpInterpolated.clear();
//		wnpInterpolated.clear();
//		unpXinterpolated.clear();
//		unpYinterpolated.clear();
//		unpZinterpolated.clear();
//
//		dphidx.clear(); dphidy.clear(); dphidz.clear();
//		d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/
//
//		// interpolate v in x
//	nn = CreateInterpolateXSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpInterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(vnpInterpolated, nn, aa, bb, cc);
//
//	// interpolate w in x
//	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpInterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(wnpInterpolated, nn, aa, bb, cc);
//
//	// interpolate u in x, y, z
//
//	nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpXinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(unpXinterpolated, nn, aa, bb, cc);
//
//	nn = CreateInterpolateYSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpYinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(unpYinterpolated, nn, aa, bb, cc);
//
//	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpZinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(unpZinterpolated, nn, aa, bb, cc);
//
//	// u mom eq
//
//	// calc duudx
//	nn = CreatedphidxUConvSystemTD(unpXinterpolated, aa, bb, cc, dphidx);
//
//	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);
//
//	// calc duvdy
//	nn = CreatedphidyUConvSystemTD(unpYinterpolated, vnpInterpolated, aa, bb, cc, dphidy);
//
//	LSSolver::SolveTridiagonalDestructive(dphidy, nn, aa, bb, cc);
//
//	// calc duwdz
//	nn = CreatedphidzUConvSystemTD(unpZinterpolated, wnpInterpolated, aa, bb, cc, dphidz);
//
//	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);
//
//	// calc d2udx2
//	nn = Created2phidx2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, d2phidx2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);
//
//	// calc d2udy2
//	nn = Created2phidy2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, d2phidy2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidy2, nn, aa, bb, cc);
//
//	// calc d2udz2
//	nn = Created2phidz2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, d2phidz2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);
//
//	// get dudt
//
//	for (int j = 1; j < cNYU - 1; j++)
//	{
//		for (int k = 1; k < cNZU - 1; k++)
//		{
//			for (int i = 1; i < cNXU - 1; i++)
//			{
//				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
//				{
//					dphidt[i + (k*cNXU) + (j* cNXU* cNZU)] = 0.0;
//				}
//				else
//				{
//					dphidt[i + (k*cNXU) + (j* cNXU* cNZU)] =
//						-1.0*(dphidx[i + (k* cNXU) + (j* cNXU* cNZU)] +
//							dphidy[i + (k* cNXU) + (j* cNXU* cNZU)] +
//							dphidz[i + (k* cNXU) + (j* cNXU* cNZU)])
//						+ (cmu / crho)*(d2phidx2[i + (k* cNXU) + (j* cNXU* cNZU)]
//							+ d2phidy2[i + (k* cNXU) + (j* cNXU* cNZU)] +
//							d2phidz2[i + (k* cNXU) + (j* cNXU* cNZU)]);
//				}
//			}
//		}
//	}
//}
//
//void SquareCylinder3DStaggeredField::CalcdphidtVmomDiv(field3D& dphidt,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc)
//{
//	int nn = 0;
//
//	static Field3DJIK unpInterpolated;
//	static Field3DJIK wnpInterpolated;
//	static Field3DIKJ vnpXinterpolated;
//	static Field3DJIK vnpYinterpolated;
//	static Field3DKJI vnpZinterpolated;
//
//	static field3D dphidx, dphidy, dphidz;
//	static Field3DIKJ d2phidx2;
//	static Field3DJIK d2phidy2;
//	static Field3DKJI d2phidz2;
//
//	/*unpInterpolated.clear();
//	wnpInterpolated.clear();
//	vnpXinterpolated.clear();
//	vnpYinterpolated.clear();
//	vnpZinterpolated.clear();
//
//	dphidx.clear(); dphidy.clear(); dphidz.clear();
//	d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/
//
//	// interpolate u in y
//	nn = CreateInterpolateYSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpInterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(unpInterpolated, nn, aa, bb, cc);
//
//	// interpolate w in y
//
//	nn = CreateInterpolateYSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpInterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(wnpInterpolated, nn, aa, bb, cc);
//
//	// interpolate v in x, y, z
//
//	nn = CreateInterpolateXSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpXinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(vnpXinterpolated, nn, aa, bb, cc);
//
//	nn = CreateInterpolateYSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpYinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(vnpYinterpolated, nn, aa, bb, cc);
//
//	nn = CreateInterpolateZSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpZinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(vnpZinterpolated, nn, aa, bb, cc);
//
//	// v mom eq
//
//	// calc duvdx
//	nn = CreatedphidxVConvSystemTD(vnpXinterpolated, unpInterpolated, aa, bb, cc, dphidx);
//
//	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);
//
//	// calc dvvdy
//	nn = CreatedphidyVConvSystemTD(vnpYinterpolated, aa, bb, cc, dphidy);
//
//	LSSolver::SolveTridiagonalDestructive(dphidy, nn, aa, bb, cc);
//
//	// calc dvwdz
//	nn = CreatedphidzVConvSystemTD(vnpZinterpolated, wnpInterpolated, aa, bb, cc, dphidz);
//
//	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);
//
//	// calc d2vdx2
//	nn = Created2phidx2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, d2phidx2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);
//
//	// calc d2vdy2
//	nn = Created2phidy2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, d2phidy2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidy2, nn, aa, bb, cc);
//
//	// calc d2vdz2
//	nn = Created2phidz2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, d2phidz2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);
//
//	// get dvdt
//
//#pragma omp parallel for
//	for (int j = 1; j < cNYV - 1; j++)
//	{
//		for (int k = 1; k < cNZV - 1; k++)
//		{
//			for (int i = 1; i < cNXV - 1; i++)
//			{
//				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd  &&
//					j >= cjCyInit - 1 && j < cjCyEnd)
//				{
//					dphidt[i + (k*cNXV) + (j* cNXV* cNZV)] = 0.0;
//				}
//				else
//				{
//					dphidt[i + (k*cNXV) + (j* cNXV* cNZV)] =
//						-1.0*(dphidx[i + (k* cNXV) + (j* cNXV* cNZV)] +
//							dphidy[i + (k* cNXV) + (j* cNXV* cNZV)] +
//							dphidz[i + (k* cNXV) + (j* cNXV* cNZV)])
//						+ (cmu / crho)*(d2phidx2[i + (k* cNXV) + (j* cNXV* cNZV)]
//							+ d2phidy2[i + (k* cNXV) + (j* cNXV* cNZV)] +
//							d2phidz2[i + (k* cNXV) + (j* cNXV* cNZV)]);
//				}
//			}
//		}
//	}
//}
//
//void SquareCylinder3DStaggeredField::CalcdphidtWmomDiv(field3D& dphidt,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc)
//{
//	int nn = 0;
//
//	static field3D unpInterpolated;
//	static field3D vnpInterpolated;
//	static field3D wnpXinterpolated;
//	static field3D wnpYinterpolated;
//	static field3D wnpZinterpolated;
//
//	static field3D dphidx, dphidy, dphidz;
//	static field3D d2phidx2, d2phidy2, d2phidz2;
//
//	/*unpInterpolated.clear();
//	vnpInterpolated.clear();
//	wnpXinterpolated.clear();
//	wnpYinterpolated.clear();
//	wnpZinterpolated.clear();
//
//	dphidx.clear(); dphidy.clear(); dphidz.clear();
//	d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/
//
//	// interpolate u in z
//	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpInterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(unpInterpolated, nn, aa, bb, cc);
//
//	// interpolate v in z
//
//	nn = CreateInterpolateZSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpInterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(vnpInterpolated, nn, aa, bb, cc);
//
//	// interpolate w in x, y, z
//
//	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpXinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(wnpXinterpolated, nn, aa, bb, cc);
//
//	nn = CreateInterpolateYSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpYinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(wnpYinterpolated, nn, aa, bb, cc);
//
//	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpZinterpolated);
//
//	LSSolver::SolveTridiagonalDestructive(wnpZinterpolated, nn, aa, bb, cc);
//
//	// w mom eq
//
//	// calc duwdx
//	nn = CreatedphidxWConvSystemTD(wnpXinterpolated, unpInterpolated, aa, bb, cc, dphidx);
//
//	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);
//
//	// calc dvwdy
//	nn = CreatedphidyWConvSystemTD(wnpYinterpolated, vnpInterpolated, aa, bb, cc, dphidy);
//
//	LSSolver::SolveTridiagonalDestructive(dphidy, nn, aa, bb, cc);
//
//	// calc dwwdz
//	nn = CreatedphidzWConvSystemTD(wnpZinterpolated, aa, bb, cc, dphidz);
//
//	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);
//
//	// calc d2udx2
//	nn = Created2phidx2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, d2phidx2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);
//
//	// calc d2udy2
//	nn = Created2phidy2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, d2phidy2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidy2, nn, aa, bb, cc);
//
//	// calc d2udz2
//	nn = Created2phidz2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, d2phidz2);
//
//	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);
//
//	// get dudt
//
//	for (int j = 0; j < cNYW; j++)
//	{
//		for (int k = 0; k < cNZW; k++)
//		{
//			for (int i = 0; i < cNXW; i++)
//			{
//				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
//					j >= cjCyInit && j < cjCyEnd)
//				{
//					dphidt[i + (k*cNXW) + (j* cNXW* cNZW)] = 0.0;
//				}
//				else
//				{
//					dphidt[i + (k*cNXW) + (j* cNXW* cNZW)] =
//						-1.0*(dphidx[i + (k* cNXW) + (j* cNXW* cNZW)] +
//							dphidy[i + (k* cNXW) + (j* cNXW* cNZW)] +
//							dphidz[i + (k* cNXW) + (j* cNXW* cNZW)])
//						+ (cmu / crho)*(d2phidx2[i + (k* cNXW) + (j* cNXW* cNZW)]
//							+ d2phidy2[i + (k* cNXW) + (j* cNXW* cNZW)] +
//							d2phidz2[i + (k* cNXW) + (j* cNXW* cNZW)]);
//				}
//			}
//		}
//	}
//}


void SquareCylinder3DStaggeredField::CalcdphidtUmomAdv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	/*static Field3DJIK vnpc;
	static Field3DKJI wnpc;*/

	static Field3DIKJ vnph;
	static Field3DIKJ wnph;

	static Field3DIKJ dphidx;
	static Field3DJIK dphidy;
	static Field3DKJI dphidz;

	static Field3DIKJ d2phidx2;
	static Field3DJIK d2phidy2;
	static Field3DKJI d2phidz2;
	/*
		vnpInterpolated.clear();
		wnpInterpolated.clear();
		unpXinterpolated.clear();
		unpYinterpolated.clear();
		unpZinterpolated.clear();

		dphidx.clear(); dphidy.clear(); dphidz.clear();
		d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();
	*/

	//// interpolate v in y
	//nn = CreateInterpolateYSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpc);

	//LSSolver::SolveTridiagonalDestructive(vnpc, nn, aa, bb, cc);

	//// interpolate w in z
	//nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpc);

	//LSSolver::SolveTridiagonalDestructive(wnpc, nn, aa, bb, cc);

	// interpolate vc in x
	nn = CreateInterpolateXPreIntSystemTD(cVnpc, cNXV, cNYV - 1, cNZV, 1, aa, bb, cc, vnph);

	LSSolver::SolveTridiagonalDestructive(vnph, nn, aa, bb, cc);

	// interpolate wc in x
	nn = CreateInterpolateXPreIntSystemTD(cWnpc, cNXW, cNYW, cNZW - 1, 2, aa, bb, cc, wnph);

	LSSolver::SolveTridiagonalDestructive(wnph, nn, aa, bb, cc);

	/*SaveIKCutToCSV("vnp", cVnp, cNXV, (cNYV), cNZV, cNYP / 2);	
	
	SaveIKCutToCSV("vnph", vnph, cNXV - 1, (cNYV - 1), cNZV, cNYP / 2);
	
	SaveIKCutToCSV("wnp", cWnp, cNXW, cNYW, (cNZW), cNYP / 2);
		
	SaveIKCutToCSV("wnph", wnph, cNXW - 1, cNYW, (cNZW - 1), cNYP / 2);*/
	

	// u mom eq

	// calc dudx
	//nn = CreatedphidxSystemTDSUCS(cUnp, cNXU, cNYU, cNZU, 0, cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDOUCS3(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, dphidx);
	nn = CreatedphidxSystemTDUPWIND(cUnp, cNXU, cNYU, cNZU, 0, cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);

	// calc dudy
	//nn = CreatedphidySystemTDSUCS(cUnp, cNXU, cNYU, cNZU, 0, vnph, (cNXV - 1), (cNYV - 1), cNZV, 1, aa, bb, cc, dphidy);
	//nn = CreatedphidySystemTDOUCS3(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, dphidy);
	nn = CreatedphidySystemTDUPWIND(cUnp, cNXU, cNYU, cNZU, 0, vnph, (cNXV - 1), (cNYV - 1), cNZV, 1, aa, bb, cc, dphidy);

	LSSolver::SolveTridiagonalDestructive(dphidy, nn, aa, bb, cc);


	// calc dudz
	//nn = CreatedphidzSystemTDSUCS(cUnp, cNXU, cNYU, cNZU, 0, wnph, (cNXW - 1), cNYW, cNZW - 1, 2, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDOUCS3(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, dphidz);
	nn = CreatedphidzSystemTDUPWIND(cUnp, cNXU, cNYU, cNZU, 0, wnph, (cNXW - 1), cNYW, cNZW - 1, 2, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	/*SaveIKCutToCSV("unp", cUnp, cNXU, cNYU, cNZU, cNYU / 2);

	SaveIKCutToCSV("dudx", dphidx, cNXU, cNYU, (cNZU), cNYU / 2);
	SaveIKCutToCSV("dudy", dphidy, cNXU, cNYU, (cNZU), cNYU / 2);
	SaveIKCutToCSV("dudz", dphidz, cNXU, cNYU, (cNZU), cNYU / 2);*/

	// calc d2udx2
	nn = Created2phidx2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);

	// calc d2udy2
	nn = Created2phidy2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, d2phidy2);

	LSSolver::SolveTridiagonalDestructive(d2phidy2, nn, aa, bb, cc);

	// calc d2udz2
	nn = Created2phidz2SystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	// get dudt
#pragma omp parallel for
	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(
							cUnp(i, j, k) * dphidx(i, j, k) +

							vnph(i, j-1, k) * dphidy(i, j, k) +

							wnph(i, j, k-1) * dphidz(i, j, k))

						+ (cmu / crho)*(d2phidx2(i, j, k)
							+ d2phidy2(i, j, k) +
							d2phidz2(i, j, k));
				}
			}
		}
	}
}

void SquareCylinder3DStaggeredField::CalcdphidtVmomAdv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	/*static Field3DIKJ unpc;
	static Field3DKJI wnpc;*/
	static Field3DJIK unpv;
	static Field3DJIK wnpv;

	static Field3DIKJ dphidx;
	static Field3DJIK dphidy;
	static Field3DKJI dphidz;
	static Field3DIKJ d2phidx2;
	static Field3DJIK d2phidy2;
	static Field3DKJI d2phidz2;

	/*unpInterpolated.clear();
	wnpInterpolated.clear();
	vnpXinterpolated.clear();
	vnpYinterpolated.clear();
	vnpZinterpolated.clear();

	dphidx.clear(); dphidy.clear(); dphidz.clear();
	d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/

	//// interpolate u in x
	//nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpc);

	//LSSolver::SolveTridiagonalDestructive(unpc, nn, aa, bb, cc);

	//// interpolate w in z

	//nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, wnpc);

	//LSSolver::SolveTridiagonalDestructive(wnpc, nn, aa, bb, cc);

	// interpolate u in y

	nn = CreateInterpolateYPreIntSystemTD(cUnpc, (cNXU - 1), cNYU, cNZU, 0, aa, bb, cc, unpv);

	LSSolver::SolveTridiagonalDestructive(unpv, nn, aa, bb, cc);

	// interpolate w in y

	nn = CreateInterpolateYPreIntSystemTD(cWnpc, cNXW, cNYW, (cNZW - 1), 2, aa, bb, cc, wnpv);

	LSSolver::SolveTridiagonalDestructive(wnpv, nn, aa, bb, cc);

	//SaveIKCutToCSV("unp", cUnp, cNXU, (cNYU), cNZU, cNYP / 2);

	//SaveIKCutToCSV("unpv", unpv, cNXU - 1, cNYU-1, cNZU, cNYP / 2);
	

	// v mom eq

	// calc dvdx
	//nn = CreatedphidxSystemTDSUCS(cVnp, cNXV, cNYV, cNZV, 1, unpv, cNXU - 1, (cNYU - 1), cNZU, 0, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDOUCS3(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, dphidx);
	nn = CreatedphidxSystemTDUPWIND(cVnp, cNXV, cNYV, cNZV, 1, unpv, cNXU - 1, (cNYU - 1), cNZU, 0, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);

	// calc dvdy
	//nn = CreatedphidySystemTDSUCS(cVnp, cNXV, cNYV, cNZV, 1, cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, dphidy);
	//nn = CreatedphidySystemTDOUCS3(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, dphidy);
	nn = CreatedphidySystemTDUPWIND(cVnp, cNXV, cNYV, cNZV, 1, cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, dphidy);

	LSSolver::SolveTridiagonalDestructive(dphidy, nn, aa, bb, cc);

	// calc dvdz
	//nn = CreatedphidzSystemTDSUCS(cVnp, cNXV, cNYV, cNZV, 1, wnpv, cNXW, cNYW - 1, cNZW - 1, 2, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDOUCS3(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, dphidz);
	nn = CreatedphidzSystemTDUPWIND(cVnp, cNXV, cNYV, cNZV, 1, wnpv, cNXW, cNYW - 1, cNZW - 1, 2, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	// calc d2vdx2
	nn = Created2phidx2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);

	// calc d2vdy2
	nn = Created2phidy2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, d2phidy2);

	LSSolver::SolveTridiagonalDestructive(d2phidy2, nn, aa, bb, cc);

	// calc d2vdz2
	nn = Created2phidz2SystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	// get dvdt

#pragma omp parallel for
	for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			for (int i = 1; i < cNXV - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd  &&
					j >= cjCyInit - 1 && j < cjCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(
							unpv(i-1, j, k) * dphidx(i, j, k) +
							cVnp(i, j, k) * dphidy(i, j, k) +
							wnpv(i, j, k-1) * dphidz(i, j, k))

						+ (cmu / crho)*(d2phidx2(i, j, k)
							+ d2phidy2(i, j, k) +
							d2phidz2(i, j, k));
				}
			}
		}
	}
}

void SquareCylinder3DStaggeredField::CalcdphidtWmomAdv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	//static Field3DIKJ unpc;
	//static Field3DJIK vnpc;
	static Field3DKJI unpt;
	static Field3DKJI vnpt;

	static Field3DIKJ dphidx;
	static Field3DJIK dphidy;
	static Field3DKJI dphidz;
	static Field3DIKJ d2phidx2;
	static Field3DJIK d2phidy2;
	static Field3DKJI d2phidz2;

	/*unpInterpolated.clear();
	vnpInterpolated.clear();
	wnpXinterpolated.clear();
	wnpYinterpolated.clear();
	wnpZinterpolated.clear();

	dphidx.clear(); dphidy.clear(); dphidz.clear();
	d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/

	//// interpolate u in x
	//nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNYU, cNZU, 0, aa, bb, cc, unpc);

	//LSSolver::SolveTridiagonalDestructive(unpc, nn, aa, bb, cc);

	//// interpolate v in y

	//nn = CreateInterpolateYSystemTD(cVnp, cNXV, cNYV, cNZV, 1, aa, bb, cc, vnpc);

	//LSSolver::SolveTridiagonalDestructive(vnpc, nn, aa, bb, cc);

	// interpolate u in z

	nn = CreateInterpolateZPreIntSystemTD(cUnpc, (cNXU - 1), cNYU, cNZU, 0, aa, bb, cc, unpt);

	LSSolver::SolveTridiagonalDestructive(unpt, nn, aa, bb, cc);

	// interpolate v in z

	nn = CreateInterpolateZPreIntSystemTD(cVnpc, cNXV, cNYV - 1, cNZV, 1, aa, bb, cc, vnpt);

	LSSolver::SolveTridiagonalDestructive(vnpt, nn, aa, bb, cc);
			
	/*SaveIKCutToCSV("vnpt", vnpt, cNXV, (cNYV - 1), cNZV - 1, cNYP / 2);
			
	SaveIKCutToCSV("unpt", unpt, cNXU-1, cNYU , cNZU - 1, cNYP / 2);*/
		
	// w mom eq

	// calc dwdx
	//nn = CreatedphidxSystemTDSUCS(cWnp, cNXW, cNYW, cNZW, 2, unpt, cNXU - 1, cNYU, cNZU - 1, 0, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDOUCS3(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, dphidx);
	nn = CreatedphidxSystemTDUPWIND(cWnp, cNXW, cNYW, cNZW, 2, unpt, cNXU - 1, cNYU, cNZU - 1, 0, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);

	// calc dwdy
	//nn = CreatedphidySystemTDSUCS(cWnp, cNXW, cNYW, cNZW, 2, vnpt, cNXV, cNYV - 1, cNZV - 1, 1, aa, bb, cc, dphidy);
	//nn = CreatedphidySystemTDOUCS3(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, dphidy);
	nn = CreatedphidySystemTDUPWIND(cWnp, cNXW, cNYW, cNZW, 2, vnpt, cNXV, cNYV - 1, cNZV - 1, 1, aa, bb, cc, dphidy);

	LSSolver::SolveTridiagonalDestructive(dphidy, nn, aa, bb, cc);

	// calc dwdz
	//nn = CreatedphidzSystemTDSUCS(cWnp, cNXW, cNYW, cNZW, 2, cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDOUCS3(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, dphidz);
	nn = CreatedphidzSystemTDUPWIND(cWnp, cNXW, cNYW, cNZW, 2, cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);

	// calc d2udy2
	nn = Created2phidy2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, d2phidy2);

	LSSolver::SolveTridiagonalDestructive(d2phidy2, nn, aa, bb, cc);

	// calc d2udz2
	nn = Created2phidz2SystemTD(cWnp, cNXW, cNYW, cNZW, 2, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	//SaveIKCutToCSV("dwdx", dphidx, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("dwdy", dphidy, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("dwdz", dphidz, cNXW, cNYW, cNZW, cNYW / 2);

	//SaveIKCutToCSV("d2wdx2", d2phidx2, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("d2wdy2", d2phidy2, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("d2wdz2", d2phidz2, cNXW, cNYW, cNZW, cNYW / 2);
	
	// get dwdt
#pragma omp parallel for
	for (int j = 1; j < cNYW - 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			for (int i = 1; i < cNXW - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					j >= cjCyInit && j < cjCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(
							unpt(i-1, j, k) * dphidx(i, j, k) +
							vnpt(i, j-1, k) * dphidy(i, j, k) +
							cWnp(i, j, k) * dphidz(i, j, k))
						+ (cmu / crho)*(d2phidx2(i, j, k)
							+ d2phidy2(i, j, k) +
							d2phidz2(i, j, k));
				}
			}
		}
	}

	/*SaveIKCutToCSV("dwdt", dphidt, cNXW, cNYW, cNZW, cNYW / 2);*/
}

void SquareCylinder3DStaggeredField::SaveIKCutToCSV(const std::string& nomeArq, const Field3D& campo, int NX, int NY, int NZ, int j)
{
	std::ofstream arq(nomeArq + ".csv");

	for (int k = 0; k < NZ; k++)
	{
		for (int i = 0; i < NX; i++)
		{
			arq << campo(i, j, k) << "\t";
		}

		arq << "\n";
	}

	arq.close();
}

void SquareCylinder3DStaggeredField::SaveIJCutToCSV(const std::string& nomeArq, const Field3D& campo, int NX, int NY, int NZ, int k)
{
	std::ofstream arq(nomeArq + ".csv");

	for (int j = 0; j < NY; j++)
	{
		for (int i = 0; i < NX; i++)
		{
			arq << campo(i, j, k) << "\t";
		}

		arq << "\n";
	}

	arq.close();
}

void SquareCylinder3DStaggeredField::SaveIKCutToCSV(const std::string& filename, int j)
{
	std::ofstream arqU(filename + "U.csv");

	for (int k = 0; k < cNZU; k++)
	{
		for (int i = 0; i < cNXU; i++)
		{
			arqU << cUnp(i, j, k) << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "V.csv");

	for (int k = 0; k < cNZV; k++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqV << cVnp(i, j, k) << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "W.csv");

	for (int k = 0; k < cNZW; k++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqW << cWnp(i, j, k) << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "P.csv");

	for (int k = 0; k < cNZP; k++)
	{
		for (int i = 0; i < cNXP; i++)
		{
			arqP << cPn(i, j, k) << "\t";
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

void SquareCylinder3DStaggeredField::SaveJKCutToCSV(const std::string& filename, const int i)
{
	std::ofstream arqU(filename + "TransU.csv");

	for (int j = 0; j < cNYU; j++)
	{
		for (int k = 0; k < cNZU; k++)
		{
			arqU << cUnp(i, j, k) << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "TransV.csv");

	for (int j = 0; j < cNYV; j++)
	{
		for (int k = 0; k < cNZV; k++)
		{
			arqV << cVnp(i, j, k) << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "TransW.csv");

	for (int j = 0; j < cNYW; j++)
	{
		for (int k = 0; k < cNZW; k++)
		{
			arqW << cWnp(i, j, k) << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "TransP.csv");

	for (int j = 0; j < cNYP; j++)
	{
		for (int k = 0; k < cNZP; k++)
		{
			arqP << cPn(i, j, k) << "\t";
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

void SquareCylinder3DStaggeredField::SaveIJCutToCSV(const std::string& filename, int k)
{
	std::ofstream arqU(filename + "TransKU.csv");

	for (int j = 0; j < cNYU; j++)
	{
		for (int i = 0; i < cNXU; i++)
		{
			arqU << cUnp(i, j, k) << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqV(filename + "TransKV.csv");

	for (int j = 0; j < cNYV; j++)
	{
		for (int i = 0; i < cNXV; i++)
		{
			arqV << cVnp(i, j, k) << "\t";
		}

		arqV << "\n";
	}

	arqV.close();

	std::ofstream arqW(filename + "TransKW.csv");

	for (int j = 0; j < cNYW; j++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqW << cWnp(i, j, k) << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "TransKP.csv");

	for (int j = 0; j < cNYP; j++)
	{
		for (int i = 0; i < cNXP; i++)
		{
			arqP << cPn(i, j, k) << "\t";
		}

		arqP << "\n";
	}

	arqP.close();
}

//int SquareCylinder3DStaggeredField::CreatedphidxUConvSystemTD(
//	const std::vector<double>& phiInterpolateX,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXU * cNYU * cNZU;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYU; j++)
//		for (int k = 0; k < cNZU; k++)
//			for (int i = 0; i < cNXU; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (i == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (i == cNXU - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((i == ciCyInit - 1 || i == ciCyInit - 2 || i == ciCyEnd - 1 || i == ciCyEnd)
//					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uipjk = phiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//
//					const double uimjk = phiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//
//					rhs(i, j, k) = ((uipjk * uipjk)
//						- (uimjk*uimjk)
//						) /
//						(ch);
//				}
//				else if (i == 1 || i == 2 || i == cNXU - 2 || i == cNXU - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uipjk = phiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//
//					const double uimjk = phiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//
//					rhs(i, j, k) = ((uipjk * uipjk)
//						- (uimjk*uimjk)
//						) /
//						(ch);
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uipjk = 0.5*(
//						cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] + cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)]
//						);
//
//					const double uimjk = 0.5*(
//						cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] + cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)]
//						);
//
//					rhs(i, j, k) = ((uipjk * uipjk)
//						- (uimjk*uimjk)
//						) /
//						(ch);
//				}*/
//				/*else if (i == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uippjk = cUnp[i + 2 + (k * cNXU) + (j * cNXU * cNZU)];
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * uijk * uijk +
//						4 * uipjk * uipjk -
//						 uippjk* uippjk);
//				}
//				else if (i == cNXU - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uimjk = cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uimmjk = cUnp[i - 2 + (k * cNXU) + (j * cNXU * cNZU)];
//
//					rhs(i, j, k) = 0.5 / ch * (3 * uijk*uijk -
//						4 * uimjk*uimjk +
//						uimmjk*uimmjk);
//				}
//				else if (i == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uimjk = cUnp[i -1 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uippjk = cUnp[i + 2 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uipppjk = cUnp[i + 3 + (k * cNXU) + (j * cNXU * cNZU)];
//
//					rhs(i, j, k) = (1.0 / ch) * ((2 * betaj2 / 3 - 1.0 / 3.0)*uimjk*uimjk -
//						(8 * betaj2 / 3 + 0.5)*uijk*uijk
//						+ (4 * betaj2 + 1.0)*uipjk*uipjk -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*uippjk +
//						(2 * betaj2 / 3.0*uipppjk*uipppjk));
//				}
//				else if (i == cNXU - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uipjk = cUnp[i + 1 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uijk = cUnp[i + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uimjk = cUnp[i - 1 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uimmjk = cUnp[i - 2 + (k * cNXU) + (j * cNXU * cNZU)];
//					const double uimmmjk = cUnp[i - 3 + (k * cNXU) + (j * cNXU * cNZU)];
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*uipjk -
//						(8 * betajm2 / 3 + 0.5)*uijk
//						+ (4 * betajm2 + 1.0)*uimjk -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*uimmjk +
//						(2 * betajm2 / 3.0*uimmmjk));
//				}*/
//				else
//				{
//
//					const double uipjk = phiInterpolateX[i + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//					const double uippjk = phiInterpolateX[i + 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//
//					const double uimjk = phiInterpolateX[i - 1 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//					const double uimmjk = phiInterpolateX[i - 2 + (k * (cNXU - 1)) + (j * (cNXU - 1) * cNZU)];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (
//						-Fimp / 6.0 * uimmjk * uimmjk
//						- Eimp / 1.0*uimjk*uimjk
//						+ Eimp / 1.0*uipjk*uipjk
//						+ Fimp / 6.0* uippjk * uippjk);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(uippjk*uippjk - uimmjk)
//						+ aII / ch * (uipjk*uipjk - uimjk * uimjk);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidyUConvSystemTD(
//	const std::vector<double>& phiInterpolateY,
//	const field3D& V,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXU * cNYU * cNZU;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYU; j++)
//		for (int k = 0; k < cNZU; k++)
//			for (int i = 0; i < cNXU; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (j == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (j == cNYU - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((j == cjCyInit || j == cjCyInit - 1 || j == cjCyEnd - 1 || j == cjCyEnd)
//					&& k >= ckCyInit && k < ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijmk = phiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
//					const double uijpk = phiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];
//
//					const double vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
//					const double vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];
//
//					rhs(i, j, k) = ((uijpk * vijpk)
//						- (uijmk * vijmk)) / ch;
//				}
//				else if (j == 1 || j == 2 || j == cNYU - 2 || j == cNYU - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijmk = phiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
//					const double uijpk = phiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];
//
//					const double vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
//					const double vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];
//
//					rhs(i, j, k) = ((uijpk * vijpk)
//						- (uijmk * vijmk)) / ch;
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijpk = 0.5*(
//						cUnp[i + (k * cNXU) + ((j) * cNXU * cNZU)] + cUnp[i + (k * cNXU) + ((j + 1) * cNXU * cNZU)]
//						);
//
//					const double uijmk = 0.5*(
//						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)]
//						);
//
//					rhs(i, j, k) = ((uijpk * V[i + (k * (cNXV-1)) + ((j) * (cNXV - 1) * cNZV)])
//						- (uijmk * V[i  + (k * (cNXV - 1)) + ((j-1) * (cNXV - 1) * cNZV)])) /ch;
//				}
//				else if (j == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijk = cUnp[i + (k*cNXU) + (j*cNXU*cNZU)];
//					const double uijpk = cUnp[i + (k*cNXU) + ((j+1)*cNXU*cNZU)];
//					const double uijppk = cUnp[i + (k*cNXU) + ((j + 2)*cNXU*cNZU)];
//
//					const double vijk = V[i + (k*(cNXV - 1)) + (j * (cNXV - 1) * cNXV)];
//					const double vijpk = V[i + (k*(cNXV - 1)) + ((j+1) * (cNXV - 1) * cNXV)];
//					const double vijppk = V[i + (k*(cNXV - 1)) + ((j+2) * (cNXV - 1) * cNXV)];
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * uijk * V[i + (k * cNXU) + (j * cNXU * cNZV)] +
//						4 * uijpk * V[i + (k * cNXU) + ((j + 1) * cNXU * cNZV)] -
//						uijppk * V[i + (k * cNXU) + ((j + 2) * cNXU * cNZV)]);
//				}
//				else if (j == cNYU - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + ((j-1) * cNXU * cNZV)] -
//						4 * cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 2) * cNXU * cNZV)] +
//						cUnp[i + (k * cNXU) + ((j - 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 3) * cNXU * cNZV)]);
//				}
//				else if (j == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 1) * cNXU * cNZV)] -
//						(8 * betaj2 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + (j * cNXU * cNZV)]
//						+ (4 * betaj2 + 1.0)*cUnp[i + (k * cNXU) + ((j + 1)* cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 1)* cNXU * cNZV)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cUnp[i + (k * cNXU) + ((j + 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 2) * cNXU * cNZV)] +
//						(2 * betaj2 / 3.0*cUnp[i + (k * cNXU) + ((j + 3) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j + 3) * cNXU * cNZV)]));
//				}
//				else if (j == cNYU - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cUnp[i + (k * cNXU) + ((j + 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j) * cNXU * cNZV)] -
//						(8 * betajm2 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * V[i + (k * cNXU) + ((j-1) * cNXU * cNZV)]
//						+ (4 * betajm2 + 1.0)*cUnp[i + (k * cNXU) + ((j - 1) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 2) * cNXU * cNZV)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cUnp[i + (k * cNXU) + ((j - 2) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 3) * cNXU * cNZV)] +
//						(2 * betajm2 / 3.0*cUnp[i + (k * cNXU) + ((j - 3) * cNXU * cNZU)] * V[i + (k * cNXU) + ((j - 4) * cNXU * cNZV)]));
//				}*/
//				else
//				{
//
//					const double uijmk = phiInterpolateY[i + (k *cNXU) + ((j - 1) * cNXU * cNZU)];
//					const double uijmmk = phiInterpolateY[i + (k *cNXU) + ((j - 2) * cNXU * cNZU)];
//
//					const double uijpk = phiInterpolateY[i + (k*cNXU) + ((j)* cNXU * cNZU)];
//					const double uijppk = phiInterpolateY[i + (k*cNXU) + ((j + 1)* cNXU * cNZU)];
//
//					const double vijpk = V[i + (k * (cNXV - 1)) + ((j) * (cNXV - 1) * cNZV)];
//					const double vijppk = V[i + (k * (cNXV - 1)) + ((j + 1) * (cNXV - 1) * cNZV)];
//
//					const double vijmk = V[i + (k * (cNXV - 1)) + ((j - 1) * (cNXV - 1) * cNZV)];
//					const double vijmmk = V[i + (k * (cNXV - 1)) + ((j - 2) * (cNXV - 1) * cNZV)];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*uijmmk * vijmmk -
//						Eimp / 1.0*uijmk * vijmk +
//						Eimp / 1.0*uijpk * vijpk +
//						Fimp / 6.0*uijppk * vijppk);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(uijppk*vijppk - uijmmk * vijmmk)
//						+ aII / ch * (uijpk*vijpk - uijmk * vijmk);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidzUConvSystemTD(
//	const std::vector<double>& phiInterpolateZ,
//	const field3D& W,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXU * cNYU * cNZU;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYU; j++)
//		for (int k = 0; k < cNZU; k++)
//			for (int i = 0; i < cNXU; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (k == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (k == cNZU - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
//					&& i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijkp = phiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
//					const double uijkm = phiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];
//
//					const double wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//					const double wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//
//					rhs(i, j, k) = ((uijkp * wijkp)
//						- (uijkm * wijkm)) /
//						(ch);
//				}
//				else if (k == 1 || k == 2 || k == cNZU - 2 || k == cNZU - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijkp = phiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
//					const double uijkm = phiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];
//
//					const double wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//					const double wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//
//					rhs(i, j, k) = ((uijkp * wijkp)
//						- (uijkm * wijkm)) /
//						(ch);
//				}
//				/*
//				else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double uijkp = 0.5*(
//						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + ((k+1) * cNXU) + ((j) * cNXU * cNZU)]
//						);
//
//					const double uijkm = 0.5*(
//						cUnp[i + (k * cNXU) + ((j)* cNXU * cNZU)] + cUnp[i + ((k-1) * cNXU) + ((j) * cNXU * cNZU)]
//						);
//
//					rhs(i, j, k) = ((uijkp * W[i + ((k) * (cNXW-1)) + ((j)* (cNXW - 1) * cNZW)])
//						- (uijkm * W[i + ((k-1) * (cNXW - 1)) + (j  * (cNXW - 1) * cNZW)])) /
//						(ch);
//				}
//				else if (k == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)] +
//						4 * cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
//						cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)]);
//				}
//				else if (k == cNZU - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)] -
//						4 * cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] +
//						cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)]);
//				}
//				else if (k == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZW)] -
//						(8 * betaj2 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)]
//						+ (4 * betaj2 + 1.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)] +
//						(2 * betaj2 / 3.0*cUnp[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZW)]));
//				}
//				else if (k == cNZU - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k) * cNXU) + ((j)* cNXU * cNZW)] -
//						(8 * betajm2 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)]
//						+ (4 * betajm2 + 1.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)] +
//						(2 * betajm2 / 3.0*cUnp[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 4) * cNXU) + ((j)* cNXU * cNZW)]));
//				}*/
//				else
//				{
//					const double uijkp = phiInterpolateZ[i + (k*cNXU) + (j*cNXU*(cNZU - 1))];
//					const double uijkpp = phiInterpolateZ[i + ((k + 1)*cNXU) + (j*cNXU*(cNZU - 1))];
//
//					const double uijkm = phiInterpolateZ[i + ((k - 1)*cNXU) + (j*cNXU*(cNZU - 1))];
//					const double uijkmm = phiInterpolateZ[i + ((k - 2)*cNXU) + (j*cNXU*(cNZU - 1))];
//
//					const double wijkp = W[i + ((k) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//					const double wijkpp = W[i + ((k + 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//
//					const double wijkm = W[i + ((k - 1) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//					const double wijkmm = W[i + ((k - 2) * (cNXW - 1)) + ((j)* (cNXW - 1) * cNZW)];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (
//						-Fimp / 6.0* uijkmm * wijkmm
//						- Eimp / 1.0* uijkm * wijkm
//						+ Eimp / 1.0*uijkp * wijkp
//						+ Fimp / 6.0*uijkpp * wijkpp);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(uijkpp*wijkpp - uijkmm * wijkmm)
//						+ aII / ch * (uijkp*wijkp - uijkm * wijkm);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidxVConvSystemTD(
//	const std::vector<double>& phiInterpolateX,
//	const field3D& U,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXV * cNYV * cNZV;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYV; j++)
//		for (int k = 0; k < cNZV; k++)
//			for (int i = 0; i < cNXV; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (i == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (i == cNXV - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
//					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vipjk = phiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//					const double vimjk = phiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//
//					const double uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
//					const double uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];
//
//					rhs(i, j, k) = ((vipjk * uipjk)
//						- (vimjk * uimjk)) /
//						(ch);
//				}
//				else if (i == 1 || i == 2 || i == cNXV - 2 || i == cNXV - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vipjk = phiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//					const double vimjk = phiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//
//					const double uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
//					const double uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];
//
//					rhs(i, j, k) = ((vipjk * uipjk)
//						- (vimjk * uimjk)) /
//						(ch);
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vipjk = 0.5*(
//						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + 1 + ((k) * cNXV) + ((j)* cNXV * cNZV)]
//						);
//
//					const double vimjk = 0.5*(
//						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i - 1 + ((k) * cNXV) + ((j)* cNXV * cNZV)]
//						);
//
//					rhs(i, j, k) = ((vipjk * U[i + (k * cNXU) + ((j)* cNXU  * cNZU)])
//						- (vimjk * U[i  - 1 + (k * cNXU) + (j  * cNXU * cNZU)])) /
//						(ch);
//				}
//				else if (i == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)] +
//						4 * cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZU)] -
//						cVnp[i + 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZU)]);
//				}
//				else if (i == cNXV - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1+ (k * cNXU) + (j * cNXU * cNZU)] -
//						4 * cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZU)] +
//						cVnp[i - 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZU)]);
//				}
//				else if (i == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZU)] -
//						(8 * betaj2 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)]
//						+ (4 * betaj2 + 1.0)*cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZU)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cVnp[i + 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZU)] +
//						(2 * betaj2 / 3.0*cVnp[i + 3 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZU)]));
//				}
//				else if (i == cNXV - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cVnp[i + 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i + (k * cNXU) + (j * cNXU * cNZU)] -
//						(8 * betajm2 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZU)]
//						+ (4 * betajm2 + 1.0)*cVnp[i - 1 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZU)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cVnp[i - 2 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZU)] +
//						(2 * betajm2 / 3.0*cVnp[i - 3 + (k * cNXV) + (j * cNXV * cNZV)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZU)]));
//				}*/
//				else
//				{
//					const double vipjk = phiInterpolateX[i + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//					const double vippjk = phiInterpolateX[i + 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//
//					const double vimjk = phiInterpolateX[i - 1 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//					const double vimmjk = phiInterpolateX[i - 2 + (k * (cNXV - 1)) + (j * (cNXV - 1) * cNZV)];
//
//					const double uipjk = U[i + (k *cNXU) + (j*cNXU*cNZU)];
//					const double uippjk = U[i + 1 + (k *cNXU) + (j*cNXU*cNZU)];
//
//					const double uimjk = U[i - 1 + (k *cNXU) + (j*cNXU*cNZU)];
//					const double uimmjk = U[i - 2 + (k *cNXU) + (j*cNXU*cNZU)];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*vimmjk * uimmjk -
//						Eimp / 1.0*vimjk * uimjk +
//						Eimp / 1.0*vipjk * uipjk +
//						Fimp / 6.0*vippjk * uippjk);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(vippjk*uippjk - vimmjk * uimmjk)
//						+ aII / ch * (vipjk*uipjk - vimjk * uimjk);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidyVConvSystemTD(
//	const std::vector<double>& phiInterpolateY,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXV * cNYV * cNZV;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYV; j++)
//		for (int k = 0; k < cNZV; k++)
//			for (int i = 0; i < cNXV; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (j == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (j == cNYV - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((j == cjCyInit - 1 || j == cjCyInit - 2 || j == cjCyEnd - 2 || j == cjCyEnd)
//					&& k >= ckCyInit && k < ckCyEnd && i >= ciCyInit && i < ciCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vijpk = phiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];
//
//					const double vijmk = phiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];
//
//					rhs(i, j, k) = ((vijpk*vijpk)
//						- (vijmk*vijmk)) /
//						(ch);
//				}
//				else if (j == 1 || j == 2 || j == cNYV - 2 || j == cNYV - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vijpk = phiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];
//
//					const double vijmk = phiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];
//
//					rhs(i, j, k) = ((vijpk*vijpk)
//						- (vijmk*vijmk)) /
//						(ch);
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vijpk = 0.5*(
//						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k)* cNXV) + ((j+1)* cNXV * cNZV)]
//						);
//
//					const double vijmk = 0.5*(
//						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k) * cNXV) + ((j-1)* cNXV * cNZV)]
//						);
//
//					rhs(i, j, k) = ((vijpk*vijpk)
//						- (vijmk*vijmk)) /
//						(ch);
//				}
//				else if (j == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] +
//						4 * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] -
//						cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)]);
//				}
//				else if (j == cNYV - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] -
//						4 * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] +
//						cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)]);
//				}
//				else if (j == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] -
//						(8 * betaj2 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)]
//						+ (4 * betaj2 + 1.0)*cVnp[i + (k * cNXV) + ((j + 1)* cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1)* cNXV * cNZV)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 2) * cNXV * cNZV)] +
//						(2 * betaj2 / 3.0*cVnp[i + (k * cNXV) + ((j + 3) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 3) * cNXV * cNZV)]));
//				}
//				else if (j == cNYV - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j + 1) * cNXV * cNZV)] -
//						(8 * betajm2 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)]
//						+ (4 * betajm2 + 1.0)*cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 1) * cNXV * cNZV)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 2) * cNXV * cNZV)] +
//						(2 * betajm2 / 3.0*cVnp[i + (k * cNXV) + ((j - 3) * cNXV * cNZV)] * cVnp[i + (k * cNXV) + ((j - 3) * cNXV * cNZV)]));
//				}*/
//				else
//				{
//					const double vijpk = phiInterpolateY[i + (k * cNXV) + (j*cNXV * cNZV)];
//					const double vijppk = phiInterpolateY[i + (k * cNXV) + ((j + 1)*cNXV * cNZV)];
//
//					const double vijmk = phiInterpolateY[i + (k * cNXV) + ((j - 1)*cNXV * cNZV)];
//					const double vijmmk = phiInterpolateY[i + (k * cNXV) + ((j - 2)*cNXV * cNZV)];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*vijmmk * vijmmk -
//						Eimp / 1.0*vijmk * vijmk +
//						Eimp / 1.0*vijpk * vijpk +
//						Fimp / 6.0*vijppk * vijppk);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(vijppk*vijppk - vijmmk * vijmmk)
//						+ aII / ch * (vijpk*vijpk - vijmk * vijmk);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidzVConvSystemTD(
//	const std::vector<double>& phiInterpolateZ,
//	const field3D& W,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXV * cNYV * cNZV;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYV; j++)
//		for (int k = 0; k < cNZV; k++)
//			for (int i = 0; i < cNXV; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (k == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (k == cNZV - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
//					&& i >= ciCyInit && i <= ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vijkp = phiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];
//
//					const double vijkm = phiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];
//
//					const double wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
//					const double wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];
//
//					rhs(i, j, k) = ((vijkp * wijkp)
//						- (vijkm * wijkm)) /
//						(ch);
//				}
//				else if (k == 1 || k == 2 || k == cNZV - 2 || k == cNZV - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vijkp = phiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];
//
//					const double vijkm = phiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];
//
//					const double wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
//					const double wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];
//
//					rhs(i, j, k) = ((vijkp * wijkp)
//						- (vijkm * wijkm)) /
//						(ch);
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double vijkp = 0.5*(
//						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k+1)* cNXV) + ((j)* cNXV * cNZV)]
//						);
//
//					const double vijkm = 0.5*(
//						cVnp[i + (k * cNXV) + ((j)* cNXV * cNZV)] + cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)]
//						);
//
//					rhs(i, j, k) = ((vijkp * W[i + ((k) * cNXW) + ((j)* cNXW * cNZW)])
//						- (vijkm * W[i + ((k - 1) * cNXW) + ((j)  * cNXW * cNZW)])) /
//						(ch);
//				}
//				else if (k == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + (k * cNXW) + (j * cNXW * cNZW)] +
//						4 * cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						cVnp[i + ((k + 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
//				}
//				else if (k == cNZV - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + ((k-1) * cNXW) + (j * cNXW * cNZW)] -
//						4 * cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
//						cVnp[i + ((k - 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]);
//				}
//				else if (k == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betaj2 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + (k * cNXW) + (j * cNXW * cNZW)]
//						+ (4 * betaj2 + 1.0)*cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cVnp[i + ((k + 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
//						(2 * betaj2 / 3.0*cVnp[i + ((k + 3) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
//				}
//				else if (k == cNZV - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cVnp[i + ((k + 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betajm2 / 3 + 0.5)*cVnp[i + (k * cNXV) + (j * cNXV * cNZV)] * W[i + ((k-1) * cNXW) + (j * cNXW * cNZW)]
//						+ (4 * betajm2 + 1.0)*cVnp[i + ((k - 1) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k-2) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cVnp[i + ((k - 2) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] +
//						(2 * betajm2 / 3.0*cVnp[i + ((k - 3) * cNXV) + ((j)* cNXV * cNZV)] * W[i + ((k - 4) * cNXW) + ((j)* cNXW * cNZW)]));
//				}*/
//				else
//				{
//					const double vijkp = phiInterpolateZ[i + ((k)* cNXV) + (j* cNXV * (cNZV - 1))];
//					const double vijkpp = phiInterpolateZ[i + ((k + 1)* cNXV) + (j* cNXV * (cNZV - 1))];
//
//					const double vijkm = phiInterpolateZ[i + ((k - 1) * cNXV) + (j* cNXV * (cNZV - 1))];
//					const double vijkmm = phiInterpolateZ[i + ((k - 2) * cNXV) + (j* cNXV * (cNZV - 1))];
//
//					const double wijkp = W[i + ((k)* cNXW) + (j*cNXW*cNZW)];
//					const double wijkpp = W[i + ((k + 1)* cNXW) + (j*cNXW*cNZW)];
//
//					const double wijkm = W[i + ((k - 1)* cNXW) + (j*cNXW*cNZW)];
//					const double wijkmm = W[i + ((k - 2)* cNXW) + (j*cNXW*cNZW)];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*vijkmm * wijkmm
//						- Eimp / 1.0* vijkm * wijkm +
//						Eimp / 1.0*vijkp * wijkp +
//						Fimp / 6.0*vijkpp * wijkpp);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(vijkpp*wijkpp - vijkmm * wijkmm)
//						+ aII / ch * (vijkp*wijkp - vijkm * wijkm);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidxWConvSystemTD(
//	const std::vector<double>& phiInterpolateX,
//	const field3D& U,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXW * cNYW * cNZW;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYW; j++)
//		for (int k = 0; k < cNZW; k++)
//			for (int i = 0; i < cNXW; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (i == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (i == cNXW - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
//					&& k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wipjk = phiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//
//					const double wimjk = phiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//
//					const double uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//					const double uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//
//					rhs(i, j, k) = ((wipjk * uipjk)
//						- (wimjk * uimjk)) /
//						(ch);
//				}
//				else if (i == 1 || i == 2 || i == cNXW - 2 || i == cNXW - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wipjk = phiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//
//					const double wimjk = phiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//
//					const double uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//					const double uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//
//					rhs(i, j, k) = ((wipjk * uipjk)
//						- (wimjk * uimjk)) /
//						(ch);
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wipjk = 0.5*(
//						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i+1 + ((k)* cNXW) + ((j)* cNXW * cNZW)]
//						);
//
//					const double wimjk = 0.5*(
//						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i-1 + ((k) * cNXW) + ((j)* cNXW * cNZW)]
//						);
//
//					rhs(i, j, k) = ((wipjk * U[i + ((k)* cNXU) + ((j)* cNXU  * (cNZU-1))])
//						- (wimjk * U[i-1 + ((k) * cNXU) + ((j)* cNXU * (cNZU - 1))])) /
//						(ch);
//				}
//				else if (i == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] +
//						4 * cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
//						cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)]);
//				}
//				else if (i == cNXW - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
//						4 * cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] +
//						cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZW)]);
//				}
//				else if (i == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
//						(8 * betaj2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)]
//						+ (4 * betaj2 + 1.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)] +
//						(2 * betaj2 / 3.0*cWnp[i + 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZW)]));
//				}
//				else if (i == cNXW - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] -
//						(8 * betajm2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i-1 + (k * cNXU) + (j * cNXU * cNZW)]
//						+ (4 * betajm2 + 1.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i -3 + (k * cNXU) + (j * cNXU * cNZW)] +
//						(2 * betajm2 / 3.0*cWnp[i - 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZW)]));
//				}*/
//				else
//				{
//					const double wipjk = phiInterpolateX[i + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//					const double wippjk = phiInterpolateX[i + 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//
//					const double wimjk = phiInterpolateX[i - 1 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//					const double wimmjk = phiInterpolateX[i - 2 + (k * (cNXW - 1)) + (j* (cNXW - 1) * cNZW)];
//
//					const double uipjk = U[i + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//					const double uippjk = U[i + 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//
//					const double uimjk = U[i - 1 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//					const double uimmjk = U[i - 2 + (k*cNXU) + (j * cNXU * (cNZU - 1))];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*wimmjk * uimmjk -
//						Eimp / 1.0*wimjk * uimjk +
//						Eimp / 1.0*wipjk * uipjk +
//						Fimp / 6.0*wippjk * uippjk);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(wippjk*uippjk - wimmjk * uimmjk)
//						+ aII / ch * (wipjk*uipjk - wimjk * uimjk);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidyWConvSystemTD(
//	const std::vector<double>& phiInterpolateY,
//	const field3D& V,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXW * cNYW * cNZW;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYW; j++)
//		for (int k = 0; k < cNZW; k++)
//			for (int i = 0; i < cNXW; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (j == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (j == cNYW - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((j == cjCyInit || j == cjCyInit - 1 || j == cjCyEnd - 1 || j == cjCyEnd)
//					&& k >= ckCyInit - 1 && k < ckCyEnd && i >= ciCyInit && i < ciCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wijpk = phiInterpolateY[i + (k * cNXW) + (j * cNXW * cNZW)];
//
//					const double wijmk = phiInterpolateY[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)];
//
//					const double vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
//					const double vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];
//
//					rhs(i, j, k) = ((wijpk * vijpk)
//						- (wijmk * vijmk)) /
//						(ch);
//				}
//				else if (j == 1 || j == 2 || j == cNYW - 2 || j == cNYW - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wijpk = phiInterpolateY[i + (k * cNXW) + (j * cNXW * cNZW)];
//
//					const double wijmk = phiInterpolateY[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)];
//
//					const double vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
//					const double vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];
//
//					rhs(i, j, k) = ((wijpk * vijpk)
//						- (wijmk * vijmk)) /
//						(ch);
//				}
//				/*
//				else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wijpk = 0.5*(
//						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k)* cNXW) + ((j+1)* cNXW * cNZW)]
//						);
//
//					const double wijmk = 0.5*(
//						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k)* cNXW) + ((j-1)* cNXW * cNZW)]
//						);
//
//					rhs(i, j, k) = ((wijpk * V[i + ((k)* cNXV) + ((j)* cNXV * (cNZV - 1))])
//						- (wijmk * V[i + ((k)* cNXV) + ((j-1)* cNXV * (cNZV-1))])) /
//						(ch);
//				}
//				else if (j == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + (j * cNXV * cNZW)] +
//						4 * cWnp[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 1) * cNXV * cNZW)] -
//						cWnp[i + (k * cNXW) + ((j + 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 2) * cNXV * cNZW)]);
//				}
//				else if (j == cNYW - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]* V[i + (k * cNXV) + ((j-1) * cNXV * cNZW)] -
//						4 * cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 2) * cNXV * cNZW)] +
//						cWnp[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 3) * cNXV * cNZW)]);
//				}
//				else if (j == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 1) * cNXV * cNZW)] -
//						(8 * betaj2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + (j * cNXV * cNZW)]
//						+ (4 * betaj2 + 1.0)*cWnp[i + (k * cNXW) + ((j + 1)* cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 1)* cNXV * cNZW)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cWnp[i + (k * cNXW) + ((j + 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 2) * cNXV * cNZW)] +
//						(2 * betaj2 / 3.0*cWnp[i + (k * cNXW) + ((j + 3) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j + 3) * cNXV * cNZW)]));
//				}
//				else if (j == cNYW - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cWnp[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j) * cNXV * cNZW)] -
//						(8 * betajm2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * V[i + (k * cNXV) + ((j-1) * cNXV * cNZW)]
//						+ (4 * betajm2 + 1.0)*cWnp[i + (k * cNXW) + ((j - 1) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 2) * cNXV * cNZW)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cWnp[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 3) * cNXV * cNZW)] +
//						(2 * betajm2 / 3.0*cWnp[i + (k * cNXW) + ((j - 3) * cNXW * cNZW)] * V[i + (k * cNXV) + ((j - 4) * cNXV * cNZW)]));
//				}
//				*/
//				else
//				{
//					const double wijpk = phiInterpolateY[i + (k * cNXW) + (j * cNXW * cNZW)];
//					const double wijppk = phiInterpolateY[i + (k * cNXW) + ((j + 1) * cNXW * cNZW)];
//
//					const double wijmk = phiInterpolateY[i + (k* cNXW) + ((j - 1)  * cNXW * cNZW)];
//					const double wijmmk = phiInterpolateY[i + (k * cNXW) + ((j - 2) * cNXW * cNZW)];
//
//					const double vijpk = V[i + (k * cNXV) + (j * cNXV * (cNZV - 1))];
//					const double vijppk = V[i + (k * cNXV) + ((j + 1) * cNXV * (cNZV - 1))];
//
//					const double vijmk = V[i + (k * cNXV) + ((j - 1) * cNXV * (cNZV - 1))];
//					const double vijmmk = V[i + (k * cNXV) + ((j - 2) * cNXV * (cNZV - 1))];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*wijmmk * vijmmk -
//						Eimp / 1.0*wijmk * vijmk +
//						Eimp / 1.0*wijpk * vijpk +
//						Fimp / 6.0*wijppk * vijppk);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(wijppk*vijppk - wijmmk * vijmmk)
//						+ aII / ch * (wijpk*vijpk - wijmk * vijmk);
//				}
//			}
//
//	return n2;
//}
//
//int SquareCylinder3DStaggeredField::CreatedphidzWConvSystemTD(
//	const std::vector<double>& phiInterpolateZ,
//	std::vector<double>& aa,
//	std::vector<double>& bb,
//	std::vector<double>& cc,
//	std::vector<double>& rhs)const
//{
//	const int n2 = cNXW * cNYW * cNZW;
//
//	aa.clear();
//	bb.clear();
//	cc.clear();
//	rhs.clear();
//
//	if (aa.size() != n2)
//		aa.resize(n2);
//
//	if (bb.size() != n2)
//		bb.resize(n2);
//
//	if (cc.size() != n2)
//		cc.resize(n2);
//
//	if (rhs.size() != n2)
//		rhs.resize(n2);
//
//	for (int j = 0; j < cNYW; j++)
//		for (int k = 0; k < cNZW; k++)
//			for (int i = 0; i < cNXW; i++)
//			{
//				const int index = rhs.GetIndex(i, j, k);
//
//				if (k == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if (k == cNZW - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.0;
//				}
//				else if ((k == ckCyInit - 1 || k == ckCyInit - 2 || k == ckCyEnd - 1 || k == ckCyEnd)
//					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wijkp = phiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];
//
//					const double wijkm = phiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];
//
//					rhs(i, j, k) = ((wijkp*wijkp)
//						- (wijkm*wijkm)) /
//						(ch);
//				}
//				else if (k == 1 || k == 2 || k == cNZW - 2 || k == cNZW - 3)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wijkp = phiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];
//
//					const double wijkm = phiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];
//
//					rhs(i, j, k) = ((wijkp*wijkp)
//						- (wijkm*wijkm)) /
//						(ch);
//				}
//				/*else if (true)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					const double wijkp = 0.5*(
//						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k+1)* cNXW) + ((j)* cNXW * cNZW)]
//						);
//
//					const double wijkm = 0.5*(
//						cWnp[i + (k * cNXW) + ((j)* cNXW * cNZW)] + cWnp[i + ((k-1)* cNXW) + ((j)* cNXW * cNZW)]
//						);
//
//					rhs(i, j, k) = ((wijkp*wijkp)
//						- (wijkm*wijkm)) /
//						(ch);
//				}
//				else if (k == 0)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] +
//						4 * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
//				}
//				else if (k == cNZW - 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] -
//						4 * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] +
//						cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)]);
//				}
//				else if (k == 1)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betaj2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]
//						+ (4 * betaj2 + 1.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
//						(2 * betaj2 / 3.0*cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
//				}
//				else if (k == cNZW - 2)
//				{
//					aa[index] = 0.0;
//					cc[index] = 0.0;
//					bb[index] = 1.0;
//
//					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betajm2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNZW) + (j * cNXW * cNZW)]
//						+ (4 * betajm2 + 1.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW *cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
//						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
//						(2 * betajm2 / 3.0*cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]));
//				}
//				*/
//				else
//				{
//					const double wijkp = phiInterpolateZ[i + (k * cNXW) + (j * cNXW + (cNZW - 1))];
//					const double wijkpp = phiInterpolateZ[i + ((k + 1) * cNXW) + (j * cNXW + (cNZW - 1))];
//
//					const double wijkm = phiInterpolateZ[i + ((k - 1) * cNXW) + (j * cNXW + (cNZW - 1))];
//					const double wijkmm = phiInterpolateZ[i + ((k - 2) * cNXW) + (j * cNXW + (cNZW - 1))];
//
//					// im
//					//aa[index] = Dimp;
//					aa[index] = alfaII;
//
//					// i
//					bb[index] = 1.0;
//
//					// ip
//					//cc[index] = Dimp;
//					cc[index] = alfaII;
//
//					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*wijkmm * wijkmm -
//						Eimp / 1.0*wijkm * wijkm +
//						Eimp / 1.0*wijkp * wijkp +
//						Fimp / 6.0*wijkpp * wijkpp);*/
//
//					rhs(i, j, k) = bII / (3 * ch)*(wijkpp*wijkpp - wijkmm * wijkmm)
//						+ aII / ch * (wijkp*wijkp - wijkm * wijkm);
//				}
//			}
//
//	return n2;
//}

int SquareCylinder3DStaggeredField::CreateInterpolateXPreIntSystemTD(
	const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const
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
		rhs.resize(n2, NXm, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}
				else if (i == NXm - 1 || i == NXm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}
				// mantendo condições de contorno ATENÇÃO TEM QUE ESTAR ANTES DO DE BAIXO (ORDEM DOS IFS IMPORTA)
				else if (uvw == 2 && (i == ciCyInit - 1 || i == ciCyEnd - 1) && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				//else if (uvw == 2 && i >= ciCyInit - 2 && i < ciCyEnd + 1 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit-1 && k < ckCyEnd-1)
				//{
				//	aa[index] = 0.0;
				//	cc[index] = 0.0;
				//	bb[index] = 1.0;

				//	rhs(i, j, k) = 0.5*(0.0/*phi(i, j, k)*/ + phi(i+1, j, k));
				//}
				else if (uvw == 2 && i >= ciCyInit - 3 && i < ciCyEnd + 2 && j >= cjCyInit && j < cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				// mantendo condições de contorno ATENÇÃO TEM QUE ESTAR ANTES DO DE BAIXO (ORDEM DOS IFS IMPORTA)
				else if (uvw == 1 && (i == ciCyInit - 1 || i == ciCyEnd - 1) && j >= cjCyInit - 1 && j < cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				//else if (uvw == 1 && i >= ciCyInit - 2 && i < ciCyEnd + 1 && j >= cjCyInit - 1 && j < cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				//{
				//	aa[index] = 0.0;
				//	cc[index] = 0.0;
				//	bb[index] = 1.0;

				//	rhs(i, j, k) = 0.5*(0.0/*phi(i, j, k)*/ + phi(i+1, j, k));
				//}
				else if (uvw == 1 && i >= ciCyInit - 3 && i < ciCyEnd + 2 && j >= cjCyInit - 1 && j < cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i+1, j, k)) / 2.0
						+ bI * (phi(i-1, j, k) + phi(i+2, j, k)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreateInterpolateYPreIntSystemTD(
	const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs) const
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
		rhs.resize(n2, NX, NYm, NZ);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (j == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j+1, k));
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j+1, k));
				}
				// MANTENDO COND. DE CONTORNO, ORDEM IMPORTA
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd - 1 && (j == cjCyInit - 1 || j == cjCyEnd - 1) && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				//else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd - 1 && j >= cjCyInit - 2 && j < cjCyEnd + 1 && k >= ckCyInit && k < ckCyEnd)
				//{
				//	aa[index] = 0.0;
				//	cc[index] = 0.0;
				//	bb[index] = 1.0;

				//	rhs(i, j, k) = 0.5*(0.0/*phi(i, j, k)*/ + phi(i, j+1, k));
				//}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd - 1 && j >= cjCyInit - 3 && j < cjCyEnd + 2 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j + 1, k));
				}
				// MANTENDO COND. DE CONTORNO, ORDEM IMPORTA
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && (j == cjCyInit - 1 || j == cjCyEnd - 1) && k >= ckCyInit - 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				//else if (uvw == 2 && i >= ciCyInit  && i < ciCyEnd && j >= cjCyInit - 2 && j < cjCyEnd + 1 && k >= ckCyInit-1 && k < ckCyEnd-1)
				//{
				//	aa[index] = 0.0;
				//	cc[index] = 0.0;
				//	bb[index] = 1.0;

				//	rhs(i, j, k) = 0.5*(/*phi(i, j, k)*/0.0 + phi(i, j+1, k));
				//}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 3 && j < cjCyEnd + 2 && k >= ckCyInit - 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j + 1, k));
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i, j+1, k)) / 2.0
						+ bI * (phi(i, j-1, k) + phi(i, j+2, k)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreateInterpolateZPreIntSystemTD(
	const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * NY * NZm;

	/*aa.clear();
	bb.clear();
	cc.clear();
	rhs.clear();*/

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.size() != n2)
		rhs.resize(n2, NX, NY, NZm);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}
				// MANTENDO COND. CONTORNO ORDEM IMPORTA
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd - 1
					&& (k == ckCyInit - 1 || k == ckCyEnd - 1))
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				//else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd - 1 
				//	&& k >= ckCyInit-2 && k < ckCyEnd + 1)
				//{
				//	aa[index] = 0.0;
				//	cc[index] = 0.0;
				//	bb[index] = 1.0;

				//	rhs(i, j, k) = 0.5*(0.0/*phi(i, j, k)*/ + phi(i, j, k+1));
				//}
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd - 1
					&& k >= ckCyInit - 3 && k < ckCyEnd + 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				// MANTENDO COND. CONTORNO ORDEM IMPORTA
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd-1 && j >= cjCyInit && j < cjCyEnd
					&& (k == ckCyInit - 1 || k == ckCyEnd - 1))
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				//else if (uvw == 0 && i >= ciCyInit-1 && i < ciCyEnd -1 && j >= cjCyInit  && j < cjCyEnd
				//	&& k >= ckCyInit - 2 && k < ckCyEnd + 1)
				//{
				//	aa[index] = 0.0;
				//	cc[index] = 0.0;
				//	bb[index] = 1.0;

				//	rhs(i, j, k) = 0.5*(0.0/*phi(i, j, k)*/ + phi(i, j, k+1));
				//}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd - 1 && j >= cjCyInit && j < cjCyEnd
					&& k >= ckCyInit - 3 && k < ckCyEnd + 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i, j, k+1)) / 2.0
						+ bI * (phi(i, j, k-1) + phi(i, j, k+2)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreateInterpolateXSystemTD(
	const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * NY * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.cField.clear();

	rhs.cNX = NXm;
	rhs.cNY = NY;
	rhs.cNZ = NZ;

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.cField.size() != n2)
		rhs.resize(n2, NXm, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
						rhs(i, j, k) = 0.0;// 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
						rhs(i, j, k) = 0.0;// 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}				
				else if (uvw == 0 && (i >= ciCyInit - 3/* || i == ciCyEnd - 1 || i == ciCyEnd ||*/ && i < ciCyEnd+1)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}
				else if (uvw == 1 && (i >= ciCyInit - 1 && i <= ciCyEnd - 1)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0; //0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 1 && (i == ciCyInit - 2 || i == ciCyEnd || i == ciCyInit - 3 || i == ciCyEnd + 1)
					&& k >= ckCyInit && k < ckCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}
				else if (uvw == 2 && (i >= ciCyInit - 1 && i <= ciCyEnd - 1)
					&& k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 2 && (i == ciCyInit - 2 || i == ciCyEnd || i == ciCyInit - 3 || i == ciCyEnd + 1)
					&& k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i+1, j, k));
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i+1, j, k)) / 2.0
						+ bI * (phi(i-1, j, k) + phi(i+2, j, k)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreateInterpolateYSystemTD(
	const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs) const
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
		rhs.resize(n2, NX, NYm, NZ);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (j == 0 || j == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2) // apenas igual a zero se for a vel no w
						rhs(i, j, k) = 0.0; // 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j+1, k));
				}
				else if (j == NYm - 1 || j == NYm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if(uvw == 2) // apenas igual a zero se for a vel no w
						rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j+1, k));
				}				
				else if (uvw == 0 && (j >= cjCyInit - 1  && j <= cjCyEnd-1)
					&& i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}
				else if (uvw == 0 && (j == cjCyInit - 2 || j == cjCyEnd || j == cjCyInit - 3 || j == cjCyEnd + 1)
					&& i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j+1, k));
				}
				else if (uvw == 1 && (j >= cjCyInit - 3 && j < cjCyEnd +1)
					&& i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i,j, k) + phi(i, j + 1, k));
				}
				else if (uvw == 2 && (j >= cjCyInit - 1 && j <= cjCyEnd - 1)
					&& i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + (k*NX) + ((j + 1)*NX*NZ)]);
				}				
				else if (uvw == 2 && (j == cjCyInit - 2 || j == cjCyEnd || j == cjCyInit - 3 || j == cjCyEnd + 1)
					&& i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j+1, k));
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i, j+1, k)) / 2.0
						+ bI * (phi(i, j-1, k) + phi(i, j+2, k)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreateInterpolateZSystemTD(
	const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
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
		rhs.resize(n2, NX, NY, NZm);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1) // só seta igual a zero se for velocidade no y
						rhs(i, j, k) = 0.0;// 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1) // só seta igual a zero se for velocidade no y
						rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}				
				else if (uvw == 0 && (k >= ckCyInit - 1 && k <= ckCyEnd-1)
					&& i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 0 && (k == ckCyInit - 2 || k == ckCyEnd || k == ckCyInit - 3 || k == ckCyEnd + 1)
					&& i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}
				else if (uvw == 1 && (k >= ckCyInit - 1 && k == ckCyEnd-1)
					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 1 && (k == ckCyInit-2 || k == ckCyEnd || k == ckCyInit - 3 || k == ckCyEnd + 1)
					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}
				else if (uvw == 2 && (k >= ckCyInit - 3  && k < ckCyEnd + 1)
					&& i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k+1));
				}
				else
				{
					// im
					aa[index] = alfaI;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaI;

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i, j, k+1)) / 2.0
						+ bI * (phi(i, j, k-1) + phi(i, j, k+2)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidxStaggeredSystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const
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
		rhs.resize(n2, NXm, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - phi(i, j, k)) / ch;
				}
				else if (i == NXm - 1 || i == NXm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - phi(i, j, k)) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi(i+2, j, k) - phi(i-1, j, k)) / (3.0*ch)
						+ aII * (phi(i+1, j, k) - phi(i, j, k)) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidyStaggeredSystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs) const
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
		rhs.resize(n2, NX, NYm, NZ);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (j == 0 || j == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j+1, k) - phi(i, j, k)) / ch;
				}
				else if (j == NYm - 1 || j == NYm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi(i, j + 2, k) - phi(i, j - 1, k)) / (3.0*ch)
						+ aII * (phi(i, j + 1, k) - phi(i, j, k)) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidzStaggeredSystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
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
		rhs.resize(n2, NX, NY, NZm);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 || k == 1/* || true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k+1) - phi(i, j, k)) / ch;
				}
				else if (k == NZm - 1 || k == NZm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi(i, j, k + 2) - phi(i, j, k - 1)) / (3.0*ch)
						+ aII * (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidxStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const

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
		rhs.resize(n2, NXm, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - phi(i, j, k)) / ch;
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi(i + 2, j, k) - phi(i - 1, j, k)) / (3.0*ch)
						+ aII * (phi(i + 1, j, k) - phi(i, j, k)) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidyStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs) const
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
		rhs.resize(n2, NX, NYm, NZ);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (j == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j +1, k) - phi(i, j, k)) / ch;
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi(i, j + 2, k) - phi(i, j - 1, k)) / (3.0*ch)
						+ aII * (phi(i, j + 1, k) - phi(i, j, k)) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidzStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
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
		rhs.resize(n2, NX, NY, NZm);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k+1) - phi(i, j, k)) / ch;
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
				else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi(i, j, k + 2) - phi(i, j, k - 1)) / (3.0*ch)
						+ aII * (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidxStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const

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
		rhs.resize(n2, NXm, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (uvw == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - phi(i, j, k)) / ch;
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k)) / ch;
				}
				else
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
				}
				/*else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi[i + 2 + (k*NX) + (j*NX*NZ)] - phi[i - 1 + (k*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + 1 + (k*NX) + (j*NX*NZ)] - phi[i + (k*NX) + (j*NX*NZ)]) / ch;
				}*/
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidyStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs) const
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
		rhs.resize(n2, NX, NYm, NZ);

	for (int j = 0; j < NYm; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (uvw == 1 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j+1,k) - phi(i, j, k)) / ch;
				}
				else if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j + 1, k)) / ch;
				}
				else if (j == NYm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k)) / ch;
				}
				else
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - phi(i, j, k)) / ch;
				}
				/*else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi[i + (k*NX) + ((j + 2)*NX*NZ)] - phi[i + (k*NX) + ((j - 1)*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + (k*NX) + ((j + 1)*NX*NZ)] - phi[i + (k*NX) + ((j)*NX*NZ)]) / ch;
				}*/
			}

	return n2;
}

int SquareCylinder3DStaggeredField::CreatedphidzStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
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
		rhs.resize(n2, NX, NY, NZm);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (uvw == 2 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j,k+1) - phi(i, j, k)) / ch;
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k + 1)) / ch;
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k)) / ch;
				}
				else
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
				/*else
				{
					// im
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfaII;

					rhs(i, j, k) = bII * (phi[i + ((k + 2)*NX) + (j*NX*NZ)] - phi[i + ((k - 1)*NX) + (j*NX*NZ)]) / (3.0*ch)
						+ aII * (phi[i + ((k + 1)*NX) + (j*NX*NZ)] - phi[i + ((k)*NX) + (j*NX*NZ)]) / ch;
				}*/
			}

	return n2;
}

int SquareCylinder3DStaggeredField::Created2phidx2SystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
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
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				/*if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi[i - 1 + (k*NX) + (j*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + 1 + (k*NX) + (j*NX*NZ)])
						/(ch*ch);
				}
				else*/
				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be + phi at i+1;
						rhs(i, j, k) = (2.0 * -phi(i+1, j, k) - 5.0*phi(i+1, j, k) +
							4.0 * phi(i+2, j, k) - phi(i+3, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i+1, j, k) +
							4.0 * phi(i+2, j, k) - phi(i+3, j, k)) / (ch*ch);
					}

				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be + phi at i-1;
						rhs(i, j, k) = (2.0 * -phi(i-1, j, k) - 5.0*phi(i-1, j, k) +
							4.0 * phi(i-2, j, k) - phi(i-3, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i-1, j, k) +
							4.0 * phi(i-2, j, k) - phi(i-3, j, k)) / (ch*ch);
					}
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						// force value to be + phi at 0;
						rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
							- phi(i, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
							+ phi(i-1, j, k)) / (ch*ch);
					}
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
					{
						rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
							+ phi(i-1, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
							+ phi(i-1, j, k)) / (ch*ch);
					}
				}
				else if (uvw == 0 && i == ciCyInit - 2 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
						+ phi(i-1, j, k)) / (ch*ch);
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i-1, j, k) +
						4.0 * phi(i-2, j, k) - phi(i-3, j, k)) / (ch*ch);
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
						+ phi(i-1, j, k)) / (ch*ch);
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i+1, j, k) +
						4.0 * phi(i+2, j, k) - phi(i+3, j, k)) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i-1, j, k)) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i-1, j, k) - 5.0*phi(i-1, j, k) +
						4.0 * phi(i-2, j, k) - phi(i-3, j, k)) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 1 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i+1, j, k) - 5.0*phi(i+1, j, k) +
						4.0 * phi(i+2, j, k) - phi(i+3, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i-1, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i-1, j, k) - 5.0*phi(i-1, j, k) +
						4.0 * phi(i-2, j, k) - phi(i-3, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i+1, j, k) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i+1, j, k) - 5.0*phi(i+1, j, k) +
						4.0 * phi(i+2, j, k) - phi(i+3, j, k)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i+2, j, k) - 2.0*phi(i, j, k) + phi(i-2, j, k))
						+ (a / (ch*ch))*(phi(i+1, j, k) - 2.0 * phi(i, j, k) + phi(i-1, j, k));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::Created2phidy2SystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DJIK& rhs)const
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
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				/*if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (j == NY - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi[i + (k*NX) + ((j-1)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + (k*NX) + ((j+1)*NX*NZ)])
						/ (ch*ch);
				}
				else */
				if (j == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2) // apenas a velocidade em z vê velocidade nula
					{
						rhs(i, j, k) = (2.0 * -phi(i, j+1, k) - 5.0*phi(i, j+1, k) +
							4.0 * phi(i, j+2, k) - phi(i, j+3, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j+1, k) +
							4.0 * phi(i, j+2, k) - phi(i, j+3, k)) / (ch*ch);
					}

				}
				else if (j == NY - 1) 
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2) // apenas a velocidade em z vê velocidade nula
					{
						rhs(i, j, k) = (2.0 * -phi(i, j-1, k) - 5.0*phi(i, j-1, k) +
							4.0 * phi(i, j-2, k) - phi(i, j-3, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j-1, k) +
							4.0 * phi(i, j-2, k) - phi(i, j-3, k)) / (ch*ch);
					}
				}
				else if (j == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2)  // apenas a velocidade em z vê velocidade nula
					{
						rhs(i, j, k) = (phi(i, j+1, k) - 2.0 *phi(i, j, k)
							- phi(i, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j+1, k) - 2.0 *phi(i, j, k)
							+ phi(i, j-1, k)) / (ch*ch);
					}
				}
				else if (j == NY - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 2)  // apenas a velocidade em z vê velocidade nula
					{
						rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
							+ phi(i, j-1, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j+1, k) - 2.0 *phi(i, j, k)
							+ phi(i, j-1, k)) / (ch*ch);
					}
				}
				else if (uvw == 0 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i, j-1, k)) / (ch*ch);
				}
				else if (uvw == 0 && j == cjCyInit && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j-1, k) - 5.0*phi(i, j-1, k) +
						4.0 * phi(i, j-1, k) - phi(i, j-3, k)) / (ch*ch);
				}
				else if (uvw == 0 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j+1, k) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 0 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j+1, k) - 5.0*phi(i, j+1, k) +
						4.0 * phi(i, j+2, k) - phi(i, j+3, k)) / (ch*ch);
				}
				else if (uvw == 1 && j == cjCyInit - 2 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - 2.0 *phi(i, j, k)
						+ phi(i, j - 1, k)) / (ch*ch);
				}
				else if (uvw == 1 && j == cjCyInit - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j - 1, k) +
						4.0 * phi(i, j -2, k) - phi(i, j -3, k)) / (ch*ch);
				}
				else if (uvw == 1 && j == cjCyEnd && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - 2.0 *phi(i, j, k)
						+ phi(i, j - 1, k)) / (ch*ch);
				}
				else if (uvw == 1 && j == cjCyEnd - 1 && k >= ckCyInit && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j + 1, k) +
						4.0 * phi(i, j + 2, k) - phi(i, j + 3, k)) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i, j - 1, k)) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyInit && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j - 1, k) - 5.0*phi(i, j - 1, k) +
						4.0 * phi(i, j -2, k) - phi(i, j -3, k)) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyEnd && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j + 1, k) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && j == cjCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd &&
					i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j + 1, k) - 5.0*phi(i, j + 1, k) +
						4.0 * phi(i, j + 2, k) - phi(i, j + 3, k)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i, j + 2, k) - 2.0*phi(i, j, k) + phi(i, j - 2, k))
						+ (a / (ch*ch))*(phi(i, j + 1, k) - 2.0 * phi(i, j, k) + phi(i, j - 1, k));
				}
			}

	return n2;
}

int SquareCylinder3DStaggeredField::Created2phidz2SystemTD(const Field3D& phi,
	const int NX, const int NY, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
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
		rhs.resize(n2, NX, NY, NZ);

	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				/*if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi[i + ((k-1)*NX) + ((j)*NX*NZ)] - 2.0 * phi[i + (k*NX) + (j*NX*NZ)] +
						phi[i + ((k+1)*NX) + ((j)*NX*NZ)])
						/ (ch*ch);
				}
				else*/
				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)  // apenas a velocidade em y vê velocidade nula
					{
						rhs(i, j, k) = (2.0 * -phi(i,j,k+1) - 5.0*phi(i, j, k + 1) +
							4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j, k + 1) +
							4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
					}

				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)   // apenas a velocidade em y vê velocidade nula
					{
						rhs(i, j, k) = (2.0 * -phi(i, j, k - 1) - 5.0*phi(i, j, k - 1) +
							4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j, k - 1) +
							4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
					}
				}
				else if (k == 1) 
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)   // apenas a velocidade em y vê velocidade nula
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
							- phi(i, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
							+ phi(i, j, k - 1)) / (ch*ch);
					}
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1)   // apenas a velocidade em y vê velocidade nula
					{
						rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
							+ phi(i, j, k - 1)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
							+ phi(i, j, k - 1)) / (ch*ch);
					}
				}
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j, k - 1) - 5.0*phi(i, j, k - 1) +
						4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j, k + 1) - 5.0*phi(i, j, k + 1) +
						4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j, k - 1) - 5.0*phi(i, j, k - 1) +
						4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j, k + 1) - 5.0*phi(i, j, k + 1) +
						4.0 * phi(i, j, k + 2) - phi(i, j, k +3)) / (ch*ch);
				}
				else if (uvw == 1 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit - 1 && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyInit - 2 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j, k - 1) +
						4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i, j, k + 1) +
						4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd && j >= cjCyInit && j < cjCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i, j, k + 2) - 2.0*phi(i, j, k) + phi(i, j, k - 2))
						+ (a / (ch*ch))*(phi(i, j, k + 1) - 2.0 * phi(i, j, k) + phi(i, j, k - 1));
				}
			}

	return n2;
}

void SquareCylinder3DStaggeredField::EnforceContinuityVelnp(double dt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<double>& val,
	Field3DIKJ& rhs,
	bool show,
	LSSolver& lssolver)
{
	int nn = 0;

	if (true)
	{
		nn = CreatedphidxStaggeredSystemTDSimplified(cUnp, cNXU, cNYU, cNZU, 0,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		nn = CreatedphidyStaggeredSystemTDSimplified(cVnp, cNXV, cNYV, cNZV, 1,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		nn = CreatedphidzStaggeredSystemTDSimplified(cWnp, cNXW, cNYW, cNZW, 2,
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

	/*SaveIKCutToCSV("dudx", cdphidx, cNXU, (cNYU), cNZU, cNYU / 2);
	SaveIKCutToCSV("dvdy", cdphidy, cNXV, (cNYV), cNZV, cNYV / 2);
	SaveIKCutToCSV("dwdz", cdphidz, cNXW, (cNYW), cNZW, cNYW / 2);

	SaveIKCutToCSV("Unp", cUnp, cNXU, (cNYU), cNZU, cNYU / 2);
	SaveIKCutToCSV("Vnp", cVnp, cNXV, (cNYV), cNZV, cNYV / 2);
	SaveIKCutToCSV("Wnp", cWnp, cNXW, (cNYW), cNZW, cNYW / 2);*/

	nn = CreateSparsePressureEquation(cdphidx, cdphidy, cdphidz,
		dt, ptr, col, val, rhs);

	if (cFirstIteration)
	{
		lssolver.PrecondtionCRS(nn, ptr, col, val);

		cFirstIteration = false;
	}

	double errOutPEqt = 0.0;
	int nitPEqt = 0;

	//SaveIKCutToCSV("rhspress", rhs, cNXP, (cNYP), cNZP, cNYP / 2);

	lssolver.SolvePreconditionedCRS(std::move(cPn.cField), nn, ptr, col, val, rhs.cField, errOutPEqt, nitPEqt);

	if(show)
		std::cout << "Pressure eq solved its " << nitPEqt << " err " << errOutPEqt << std::endl;

	//SaveIKCutToCSV("press", cPn, cNXP, (cNYP), cNZP, cNYP / 2);

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
				cPn(0, j, k) = cPn(1, j, k);
				cPn(cNXP - 1, j, k) = cPn(cNXP - 2, j, k);
			}

		nn = CreatedphidxStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);

		// set dp/dn = 0
		for (int k = 0; k < cNZP; k++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn(i, cNYP - 1, k) = cPn(i, cNYP - 2, k);
				cPn(i, 0, k) = cPn(i, 1, k);
			}

		nn = CreatedphidyStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidy);

		lssolver.SolveTridiagonalDestructive(cdphidy, nn, aa, bb, cc);

		// set dp/dn = 0
		for (int j = 0; j < cNYP; j++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn(i, j, 0) = cPn(i, j, 1);
				cPn(i, j, cNZP - 1) = cPn(i, j, cNZP - 2);
			}

		nn = CreatedphidzStaggeredSystemTD(cPn, cNXP, cNYP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}

	//SaveIKCutToCSV("dPdx", cdphidx, cNXP, (cNYP), cNZP, cNYP / 2);
	//SaveIKCutToCSV("dPdy", cdphidy, cNXP, (cNYP), cNZP, cNYP / 2);
	//SaveIKCutToCSV("dPdz", cdphidz, cNXP, (cNYP), cNZP, cNYP / 2);


	for (int j = 1; j < cNYU - 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					continue;
				}
				cUnp(i, j, k) -=
					(dt / crho)*(cdphidx(i, j, k));

			}
		}
	}

	for (int j = 1; j < cNYV - 1; j++)
	{
		for (int k = 1; k < cNZV - 1; k++)
		{
			for (int i = 1; i < cNXV - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd &&  k >= ckCyInit && k < ckCyEnd
					&& j >= cjCyInit - 1 && j < cjCyEnd)
				{
					continue;
				}
				cVnp(i, j, k) -=
					(dt / crho)*(cdphidy(i, j, k));

			}
		}
	}

	for (int j = 1; j < cNYW - 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			for (int i = 1; i < cNXW - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd
					&& j >= cjCyInit && j < cjCyEnd)
				{
					continue;
				}
				cWnp(i, j, k) -=
					(dt / crho)*(cdphidz(i, j, k));

			}
		}
	}
}

void SquareCylinder3DStaggeredField::SaveField(const std::string& baseFileName) const
{
	std::ofstream uFile(baseFileName + "U.dat", std::ios::binary);

	uFile.write(reinterpret_cast<const char*>(&cNXU), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNYU), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNZU), sizeof(int64_t));

	uFile.write(reinterpret_cast<const char*>(cUnp.data()), cUnp.size() * sizeof(double));

	uFile.close();

	std::ofstream vFile(baseFileName + "V.dat", std::ios::binary);

	vFile.write(reinterpret_cast<const char*>(&cNXV), sizeof(int64_t));
	vFile.write(reinterpret_cast<const char*>(&cNYV), sizeof(int64_t));
	vFile.write(reinterpret_cast<const char*>(&cNZV), sizeof(int64_t));

	vFile.write(reinterpret_cast<const char*>(cVnp.data()), cVnp.size() * sizeof(double));

	vFile.close();

	std::ofstream wFile(baseFileName + "W.dat", std::ios::binary);

	wFile.write(reinterpret_cast<const char*>(&cNXW), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNYW), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNZW), sizeof(int64_t));

	wFile.write(reinterpret_cast<const char*>(cWnp.data()), cWnp.size() * sizeof(double));

	wFile.close();

	std::ofstream pFile(baseFileName + "P.dat", std::ios::binary);

	pFile.write(reinterpret_cast<const char*>(&cNXP), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNYP), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNZP), sizeof(int64_t));

	pFile.write(reinterpret_cast<const char*>(cPn.data()), cPn.size() * sizeof(double));

	pFile.close();
}

void SquareCylinder3DStaggeredField::ReadField(const std::string& basefilename)
{

	std::ifstream uFile(basefilename + "U.dat", std::ios::binary);

	uFile.read(reinterpret_cast<char*>(&cNXU), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNYU), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNZU), sizeof(int64_t));

	cUn.clear();
	cUn.resize(cNXU * cNYU * cNZU, cNXU, cNYU, cNZU);

	uFile.read(reinterpret_cast<char*>(cUn.data()), cUn.size() * sizeof(double));

	uFile.close();

	cUnp = cUn;

	std::ifstream vFile(basefilename + "V.dat", std::ios::binary);

	vFile.read(reinterpret_cast<char*>(&cNXV), sizeof(int64_t));
	vFile.read(reinterpret_cast<char*>(&cNYV), sizeof(int64_t));
	vFile.read(reinterpret_cast<char*>(&cNZV), sizeof(int64_t));

	cVn.clear();
	cVn.resize(cNXV * cNYV * cNZV, cNXV, cNYV, cNZV);

	vFile.read(reinterpret_cast<char*>(cVn.data()), cVn.size() * sizeof(double));

	vFile.close();

	cVnp = cVn;

	std::ifstream wFile(basefilename + "W.dat", std::ios::binary);

	wFile.read(reinterpret_cast<char*>(&cNXW), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNYW), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNZW), sizeof(int64_t));

	cWn.clear();
	cWn.resize(cNXW * cNYW * cNZW, cNXW, cNYW, cNZW);

	wFile.read(reinterpret_cast<char*>(cWn.data()), cWn.size() * sizeof(double));

	wFile.close();

	cWnp = cWn;

	std::ifstream pFile(basefilename + "P.dat", std::ios::binary);

	pFile.read(reinterpret_cast<char*>(&cNXP), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNYP), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNZP), sizeof(int64_t));

	cPn.clear();
	cPn.resize(cNXP * cNYP * cNZP, cNXP, cNYP, cNZP);

	pFile.read(reinterpret_cast<char*>(cPn.data()), cPn.size() * sizeof(double));

	pFile.close();
}

void SquareCylinder3DStaggeredField::InterpolateUVWInRespDirections(
	std::vector<double>& aa, std::vector<double>& bb, std::vector<double>& cc,
	std::vector<double>& aa1, std::vector<double>& bb1, std::vector<double>& cc1,
	std::vector<double>& aa2, std::vector<double>& bb2, std::vector<double>& cc2)
{
	SquareCylinder3DStaggeredField& classe = *this;
		
	
	std::thread thread1([&aa, &bb, &cc, &classe]()
		{
			int nn = 0;
			// interpolate u in x
			nn = classe.CreateInterpolateXSystemTD(classe.cUnp, classe.cNXU, classe.cNYU, classe.cNZU, 0,
				aa, bb, cc, classe.cUnpc);

			LSSolver::SolveTridiagonalDestructive(classe.cUnpc, nn, aa, bb, cc);
		});
	std::thread thread2([&aa1, &bb1, &cc1, &classe]()
		{
			int nn = 0;
			// interpolate v in y
			nn = classe.CreateInterpolateYSystemTD(classe.cVnp, classe.cNXV, classe.cNYV, classe.cNZV, 1,
				aa1, bb1, cc1, classe.cVnpc);

			LSSolver::SolveTridiagonalDestructive(classe.cVnpc, nn, aa1, bb1, cc1);
		});
	
	int nn = 0;
	// interpolate w in z
	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, 2,
		aa2, bb2, cc2, cWnpc);

	LSSolver::SolveTridiagonalDestructive(cWnpc, nn, aa2, bb2, cc2);

	thread1.join();
	thread2.join();

	/*SaveIKCutToCSV("unpc", cUnpc, cUnpc.cNX, cUnpc.cNY, cUnpc.cNZ, cNYP/2);
	SaveIKCutToCSV("vnpc", cVnpc, cVnpc.cNX, cVnpc.cNY, cVnpc.cNZ, cNYP / 2);
	SaveIKCutToCSV("wnpc", cWnpc, cWnpc.cNX, cWnpc.cNY, cWnpc.cNZ, cNYP / 2);*/
	

	/*int nn = 0;
	// interpolate u in x
	nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNYU, cNZU, 0,
		aa, bb, cc, cUnpc);

	LSSolver::SolveTridiagonalDestructive(cUnpc, nn, aa, bb, cc);

	// interpolate v in y
	nn = CreateInterpolateYSystemTD(cVnp, cNXV, cNYV, cNZV, 1,
		aa1, bb1, cc1, cVnpc);

	LSSolver::SolveTridiagonalDestructive(cVnpc, nn, aa1, bb1, cc1);

	// interpolate w in z
	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNYW, cNZW, 2,
		aa2, bb2, cc2, cWnpc);

	LSSolver::SolveTridiagonalDestructive(cWnpc, nn, aa2, bb2, cc2);*/
}