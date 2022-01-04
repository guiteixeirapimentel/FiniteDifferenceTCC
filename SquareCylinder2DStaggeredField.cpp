#include "SquareCylinder2DStaggeredField.h"
#include "LSSolver.h"

#include <string>
#include <thread>

SquareCylinder2DStaggeredField::SquareCylinder2DStaggeredField(double rho,
	double mu,
	double h,
	double l,
	double L,
	double H,
	double ufreestream)
	:
	crho(rho),
	cmu(mu),
	ch(h),
	cL(L),
	cH(H),

	cl(l),

	cFirstIteration(true),

	cDistRelFrenteRect(8),

	cNXU(int(cL / ch) + 1),
	cNZU(int(cH / ch) + 2),

	cNXV(int(cL / ch) + 2),
	cNZV(int(cH / ch) + 2),

	cNXW(int(cL / ch) + 2),
	cNZW(int(cH / ch) + 1),

	cNXP(int(cL / ch) + 2),
	cNZP(int(cH / ch) + 2),

	cNXCy(cl / ch),
	cNZCy(cNXCy),

	ciCyInit(cDistRelFrenteRect* cNXCy),
	ckCyInit(cNZP / 2 - cNZCy / 2),

	ciCyEnd(cDistRelFrenteRect* cNXCy + cNXCy),
	ckCyEnd((cNZP / 2 + cNZCy / 2) + (cNZCy % 2)),

	cUfreestream(ufreestream)
{
	cUn.resize(cNXU*cNZU, cNXU, 1, cNZU, cUfreestream);
	cWn.resize(cNXW*cNZW, cNXW, 1, cNZW, 0.0);
	cPn.resize(cNXP*cNZP, cNXP, 1, cNZP, 0.0);

	cUnp = cUn;
	cWnp = cWn;
	//cPnp = cPn;

	cPhi = cPn;

	cdphidx = cPn;
	cdphidz.resize(cNXP*cNZP, cNXP, 1, cNZP, 0.0);

	cd2phidx2.resize(cNXP*cNZP, cNXP, 1, cNZP, 0.0);
	cd2phidz2.resize(cNXP*cNZP, cNXP, 1, cNZP, 0.0);

	cdUdtn = cUn;
	cdUdt1 = cUn;
	cdUdt2 = cUn;
	cdUdt3 = cUn;

	cdWdtn = cWn;
	cdWdt1 = cWn;
	cdWdt2 = cWn;
	cdWdt3 = cWn;
}

SquareCylinder2DStaggeredField::~SquareCylinder2DStaggeredField()
{}

void SquareCylinder2DStaggeredField::SetBCFlatPlate()
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

	for (int k = ckCyInit; k < ckCyEnd; k++)
		for (int i = ciCyInit - 1; i < ciCyEnd; i++)
		{
			cUn(i, 0, k) = 0.0;
		}

		
	
	for (int k = ckCyInit - 1; k < ckCyEnd; k++)
		for (int i = ciCyInit; i < ciCyEnd; i++)
		{
			cWn(i, 0, k) = 0.0;
		}

	cUnp = cUn;
	cWnp = cWn;
}

void SquareCylinder2DStaggeredField::UpdateBCFlatPlate()
{
	// extrapolate outlet (i = NX - 1)

	for (int k = 1; k < cNZU - 1; k++)
	{
		const int i = cNXU - 1;
		const int iu = cNXU - 2;

		cUnp(iu + 1, 0, k) = cUnp(iu, 0, k) +
			cWnp(i, 0, k - 1) - cWnp(i, 0, k);
	}
}

int SquareCylinder2DStaggeredField::CreatedphidxSystemTDOUCS3(const Field3D& phi,
	const int NX, const int NZ,
	int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX,1, NZ);

	for (int j = 0; j < 1; j++)
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
							+ (4 * betaj2 + 1.0)*phi(i + 1, j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2, j, k) + (2 * betaj2 / 3.0*phi(i + 3, j, k)));
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0 && false)
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*-phi(i, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k) + (2 * betajm2 / 3.0*phi(i - 3, j, k)));
					else
						rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
							+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k) + (2 * betajm2 / 3.0*phi(i - 3, j, k)));
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k) + (2 * betajm2 / 3.0*phi(i - 3, j, k)));
				}
				else if (uvw == 0 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i - 1, j, k) + phi(i - 2, j, k));
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i - 1, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i + 1, j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2, j, k) + (2 * betaj2 / 3.0*phi(i + 3, j, k)));
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i + 1, j, k) - phi(i + 2, j, k));
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i + 1, j, k) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i - 1, j, k) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i - 2, j, k) + (2 * betajm2 / 3.0*phi(i - 3, j, k)));
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i - 1, j, k) + phi(i - 2, j, k));
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i - 1, j, k) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i + 1, j, k) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i + 2, j, k) + (2 * betaj2 / 3.0*phi(i + 3, j, k)));
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i + 1, j, k) - phi(i + 2, j, k));
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
						(qm2*phi(i - 2, j, k) + qm1 * phi(i - 1, j, k) +
							+q0 * phi(i, j, k) +
							qp1 * phi(i + 1, j, k) + qp2 * phi(i + 2, j, k));
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreatedphidzSystemTDOUCS3(const Field3D& phi,
	const int NX, const int NZ,
	int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j, k + 1) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) -
						4 * phi(i, j, k - 1) + phi(i, j, k - 2));
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) +
						4 * phi(i, j, k + 1) - phi(i, j, k + 2));
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*phi(i, j, k - 1) - (8 * betaj2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betaj2 + 1.0)*phi(i, j, k + 1) - (8 * betaj2 / 3.0 + 1.0 / 6.0)*phi(i, j, k + 2) + (2 * betaj2 / 3.0*phi(i, j, k + 3)));
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*phi(i, j, k + 1) - (8 * betajm2 / 3 + 0.5)*phi(i, j, k)
						+ (4 * betajm2 + 1.0)*phi(i, j, k - 1) - (8 * betajm2 / 3.0 + 1.0 / 6.0)*phi(i, j, k - 2) + (2 * betajm2 / 3.0*phi(i, j, k - 3)));
				}
				else if (uvw == 2 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * phi(i, j, k) - 4 * phi(i, j, k - 1) + phi(i, j, k - 2));
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * phi(i, j, k) + 4 * phi(i, j, k + 1) - phi(i, j, k + 2));
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
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

int SquareCylinder2DStaggeredField::CreatedphidxSystemTDSUCS(const Field3D& phi,
	const int NX, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				double velConv = 0.0;

				if (i == 0 || i == NX - 1 || k == 0 || k == NZ - 1)
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

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));

					/*rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;*/
				}
				else if (i == NX - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));

					/*rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;*/
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}*/
				}
				else if (i == NX - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / (ch / 2.0);
					}*/
				}
				else if (uvw == 0 && i >= ciCyInit && i < ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 1 && i >= ciCyInit + 1 && i < ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit + 1 && i < ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));

					/*rhs(i, j, k) = (0.0 - phi(i - 1, j, k)) / ch;*/
				}
				else if (uvw == 0 && i == ciCyInit - 2 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}*/
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}*/
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));

					/*rhs(i, j, k) = (phi(i + 1, j, k) - 0.0) / ch;*/
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
					//if (velConv > 0.0)
					//{
					//	rhs(i, j, k) = (phi(i, j, k) - phi(i - 1, j, k)) / ch;
					//}
					//else
					//{
					//	rhs(i, j, k) = (/*phi(i + 1, j, k)*/0.0 - phi(i, j, k)) / (ch / 2);
					//}
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i-1, j, k)
						+ 48 * phi(i-1, j, k) - 36 * phi(i-2, j, k)
						+ 16 * phi(i-3, j, k) - 3 * phi(i-4, j, k));

					//rhs(i, j, k) = (0.0 - phi(i - 1, j, k)) / (ch / 2.0);
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i+1, j, k)
						- 6 * phi(i+2, j, k) + phi(i+3, j, k));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
					}*/
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i+1, j, k)
						+ 48 * phi(i+1, j, k) - 36 * phi(i+2, j, k)
						+ 16 * phi(i+3, j, k) - 3 * phi(i+4, j, k));

					//rhs(i, j, k) = (phi(i + 1, j, k) - 0.0) / (ch / 2.0);
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
							(velPosSUCScm2 * phi(i - 2, j, k)
								+ velPosSUCScm1 * phi(i - 1, j, k)
								+ velPosSUCSc0 * phi(i, j, k)
								+ velPosSUCSc1 * phi(i + 1, j, k));

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
							(velNegSUCSc2 * phi(i + 2, j, k)
								+ velNegSUCSc1 * phi(i + 1, j, k)
								+ velNegSUCSc0 * phi(i, j, k)
								+ velNegSUCScm1 * phi(i - 1, j, k)
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

int SquareCylinder2DStaggeredField::CreatedphidzSystemTDSUCS(const Field3D& phi,
	const int NX, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				double velConv;

				if (i == 0 || i == NX - 1 || k == 0 || k == NZ - 1)
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


					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k+1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));

					//rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
				}
				else if (k == NZ - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;


					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));

					//rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j, k - 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
					}*/
				}
				else if (k == NZ - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}*/
				}
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / (ch / 2.0);
					}*/
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));

					//rhs(i, j, k) = (0.0 - phi(i, j, k - 1)) / (ch / 2.0);
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));

					//rhs(i, j, k) = (phi(i, j, k + 1) - 0.0) / (ch / 2.0);
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / (ch / 2.0);
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
					}*/
				}
				else if (uvw == 2 && k == ckCyInit - 2 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - phi(i, j, k - 1)) / ch;
					}
					else
					{
						rhs(i, j, k) = (0.0 - phi(i, j, k)) / ch;
					}*/
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));

					//rhs(i, j, k) = (0.0 - phi(i, j, k - 1)) / ch;
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));

					//rhs(i, j, k) = (phi(i, j, k + 1) - 0.0) / ch;
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch*12.0)) * (-3 * phi(i, j, k - 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));

					/*if (velConv > 0.0)
					{
						rhs(i, j, k) = (phi(i, j, k) - 0.0) / ch;
					}
					else
					{
						rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
					}*/
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
							(velPosSUCScm2 * phi(i, j, k - 2)
								+ velPosSUCScm1 * phi(i, j, k - 1)
								+ velPosSUCSc0 * phi(i, j, k)
								+ velPosSUCSc1 * phi(i, j, k + 1));

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
							(velNegSUCSc2 * phi(i, j, k + 2)
								+ velNegSUCSc1 * phi(i, j, k + 1)
								+ velNegSUCSc0 * phi(i, j, k)
								+ velNegSUCScm1 * phi(i, j, k - 1)
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

int SquareCylinder2DStaggeredField::CreatedphidxSystemTDUPWIND(const Field3D& phi,
	const int NX, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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
				else if (uvw == 0 && i >= ciCyInit && i < ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit + 1 && i < ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (uvw == 0 && i == ciCyInit - 2 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i + 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * phi(i - 1, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));
				}
				else if (uvw == 1 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 1 && i == ciCyInit && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i - 1, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (uvw == 1 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (uvw == 1 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i + 1, j, k)
						+ 48 * phi(i + 1, j, k) - 36 * phi(i + 2, j, k)
						+ 16 * phi(i + 3, j, k) - 3 * phi(i + 4, j, k));
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i - 1, j, k)
						- 6 * phi(i - 2, j, k) + phi(i - 3, j, k));
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i - 1, j, k)
						+ 48 * phi(i - 1, j, k) - 36 * phi(i - 2, j, k)
						+ 16 * phi(i - 3, j, k) - 3 * phi(i - 4, j, k));
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i + 1, j, k)
						- 6 * phi(i + 2, j, k) + phi(i + 3, j, k));
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
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
								-phi(i - 1, j, k)
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

int SquareCylinder2DStaggeredField::CreatedphidzSystemTDUPWIND(const Field3D& phi,
	const int NX, const int NZ,
	int uvw,
	const Field3D& convVel,
	const int NXc, const int NZc,
	int uvwc,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (uvw == 1 && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit + 1 && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 2 && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -0.0;
				}
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));
				}
				else if (uvw == 1 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 1 && k == ckCyInit && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k - 1)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (uvw == 1 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * -phi(i, j, k + 1)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (uvw == 1 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-3 * -phi(i, j, k)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k + 1)
						- 6 * phi(i, j, k + 2) + phi(i, j, k + 3));
				}
				else if (uvw == 2 && k == ckCyInit - 2 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-3 * phi(i, j, k + 1)
						- 10 * phi(i, j, k) + 18 * phi(i, j, k - 1)
						- 6 * phi(i, j, k - 2) + phi(i, j, k - 3));
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -(1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k - 1) - 36 * phi(i, j, k - 2)
						+ 16 * phi(i, j, k - 3) - 3 * phi(i, j, k - 4));
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (1.0 / (ch * 12.0)) * (-25 * phi(i, j, k)
						+ 48 * phi(i, j, k + 1) - 36 * phi(i, j, k + 2)
						+ 16 * phi(i, j, k + 3) - 3 * phi(i, j, k + 4));
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
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
								+phi(i, j, k + 1)
								- phi(i, j, k)
								);
					}
				}
			}

	return n2;
}


int SquareCylinder2DStaggeredField::Created2phidx2SystemTD(const Field3D& phi,
	const int NX, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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

					rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i + 1, j, k) +
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

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i + 2, j, k) - 2.0*phi(i, j, k) + phi(i - 2, j, k))
						+ (a / (ch*ch))*(phi(i + 1, j, k) - 2.0 * phi(i, j, k) + phi(i - 1, j, k));
				}
			}

	return n2;
}


int SquareCylinder2DStaggeredField::Created2phidz2SystemTD(const Field3D& phi,
	const int NX, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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
						4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
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

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i, j, k + 2) - 2.0*phi(i, j, k) + phi(i, j, k - 2))
						+ (a / (ch*ch))*(phi(i, j, k + 1) - 2.0 * phi(i, j, k) + phi(i, j, k - 1));
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreateSparsePressureEquation(const Field3D& dUdx,
	const Field3D& dWdz,
	const double dt,
	std::vector<int>& ptr,
	std::vector<int>& col,
	std::vector<double>& val,
	Field3DIKJ& rhs)const
{
	const int n2 = cNXP * 1 * cNZP; // Number of points in the grid.

	ptr.clear();
	ptr.reserve(n2 + 1);
	ptr.push_back(0);

	col.clear();
	col.reserve(n2 * 7);

	val.clear();
	val.reserve(n2 * 7);

	rhs.clear();
	rhs.resize(n2, cNXP, 1, cNZP, 0.0);

	for (int j = 0; j < 1; ++j)
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
				else if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
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
					if (k == 1 || (k == ckCyEnd && i >= ciCyInit && i < ciCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkpi
					col.push_back(rhs.GetIndex(i, j, k + 1));
					if (k == cNZP - 2 || (k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);

						valijk -= 1.0;
					}


					// jkim
					col.push_back(rhs.GetIndex(i - 1, j, k));
					if (i == 1 || (i == ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}


					// jkip
					col.push_back(rhs.GetIndex(i + 1, j, k));
					if (i == cNXP - 2 || (i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd))
					{
						val.push_back(0.0);
					}
					else
					{
						val.push_back(1.0);
						valijk -= 1.0;
					}

					//// jmki
					//col.push_back(rhs.GetIndex(i, j - 1, k));
					//if (j == 1 || (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					//{
					//	val.push_back(0.0);
					//}
					//else
					//{
					//	val.push_back(1.0);
					//	valijk -= 1.0;
					//}


					//// jpki
					//col.push_back(rhs.GetIndex(i, j + 1, k));
					//if (j == cNYP - 2 || (j == cjCyInit - 1 && i >= ciCyInit && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd))
					//{
					//	val.push_back(0.0);
					//}
					//else
					//{
					//	val.push_back(1.0);
					//	valijk -= 1.0;
					//}

					// jki
					col.push_back(rhs.GetIndex(i, j, k));
					val.push_back(valijk);

					rhs(i, j, k) = ch * ch*(crho / dt)*(
						dUdx(i - 1, j, k) +
						/*dVdy(i, j - 1, k) +*/
						dWdz(i, j, k - 1)
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

void SquareCylinder2DStaggeredField::CalcdphidtUmomDiv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	static Field3DIKJ wnpInterpolated;
	static Field3DIKJ unpXinterpolated;
	static Field3DKJI unpZinterpolated;

	static Field3DIKJ dphidx;
	static Field3DKJI dphidz;
	static Field3DIKJ d2phidx2; 
	static Field3DKJI d2phidz2;

	/*
		vnpInterpolated.clear();
		wnpInterpolated.clear();
		unpXinterpolated.clear();
		unpYinterpolated.clear();
		unpZinterpolated.clear();

		dphidx.clear(); dphidy.clear(); dphidz.clear();
		d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/
	
	// interpolate w in x
	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, wnpInterpolated);

	LSSolver::SolveTridiagonalDestructive(wnpInterpolated, nn, aa, bb, cc);

	// interpolate u in x, y, z

	nn = CreateInterpolateXSystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, unpXinterpolated);

	LSSolver::SolveTridiagonalDestructive(unpXinterpolated, nn, aa, bb, cc);
	
	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, unpZinterpolated);

	LSSolver::SolveTridiagonalDestructive(unpZinterpolated, nn, aa, bb, cc);

	// u mom eq

	// calc duudx
	nn = CreatedphidxUConvSystemTD(unpXinterpolated, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);


	// calc duwdz
	nn = CreatedphidzUConvSystemTD(unpZinterpolated, wnpInterpolated, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);
	
	// calc d2udz2
	nn = Created2phidz2SystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 0; j < 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(dphidx(i, j, k) +
							dphidz(i, j, k))
						+ (cmu / crho)*(d2phidx2(i, j, k) +
							d2phidz2(i, j, k));
				}
			}
		}
	}
}

void SquareCylinder2DStaggeredField::CalcdphidtWmomDiv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	static Field3DKJI unpInterpolated;
	static Field3DIKJ wnpXinterpolated;
	static Field3DKJI wnpZinterpolated;
		   
	static Field3DIKJ dphidx;
	static Field3DKJI dphidz;
	static Field3DIKJ d2phidx2;
	static Field3DKJI d2phidz2;

	/*unpInterpolated.clear();
	vnpInterpolated.clear();
	wnpXinterpolated.clear();
	wnpYinterpolated.clear();
	wnpZinterpolated.clear();

	dphidx.clear(); dphidy.clear(); dphidz.clear();
	d2phidx2.clear(); d2phidy2.clear(); d2phidz2.clear();*/

	// interpolate u in z
	nn = CreateInterpolateZSystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, unpInterpolated);

	LSSolver::SolveTridiagonalDestructive(unpInterpolated, nn, aa, bb, cc);

	// interpolate w in x, y, z

	nn = CreateInterpolateXSystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, wnpXinterpolated);

	LSSolver::SolveTridiagonalDestructive(wnpXinterpolated, nn, aa, bb, cc);


	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, wnpZinterpolated);

	LSSolver::SolveTridiagonalDestructive(wnpZinterpolated, nn, aa, bb, cc);

	// w mom eq

	// calc duwdx
	nn = CreatedphidxWConvSystemTD(wnpXinterpolated, unpInterpolated, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);
	
	// calc dwwdz
	nn = CreatedphidzWConvSystemTD(wnpZinterpolated, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);
		
	// calc d2udz2
	nn = Created2phidz2SystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 0; j < 1; j++)
	{
		for (int k = 0; k < cNZW; k++)
		{
			for (int i = 0; i < cNXW; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(dphidx(i, j, k) +
							dphidz(i, j, k))
						+ (cmu / crho)*(d2phidx2(i, j, k) +
							d2phidz2(i, j, k));
				}
			}
		}
	}
}


void SquareCylinder2DStaggeredField::CalcdphidtUmomAdv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	/*
	static Field3DKJI wnpc;*/
	static Field3DIKJ wnph;

	static Field3DIKJ dphidx;
	static Field3DKJI dphidz;

	static Field3DIKJ d2phidx2;
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
	
	// interpolate wc in x
	nn = CreateInterpolateXPreIntSystemTD(cWnpc, cNXW, cNZW - 1, 2, aa, bb, cc, wnph);

	LSSolver::SolveTridiagonalDestructive(wnph, nn, aa, bb, cc);

	/*SaveIKCutToCSV("vnp", cVnp, cNXV, (cNYV), cNZV, cNYP / 2);

	SaveIKCutToCSV("vnph", vnph, cNXV - 1, (cNYV - 1), cNZV, cNYP / 2);

	SaveIKCutToCSV("wnp", cWnp, cNXW, cNYW, (cNZW), cNYP / 2);

	SaveIKCutToCSV("wnph", wnph, cNXW - 1, cNYW, (cNZW - 1), cNYP / 2);*/


	// u mom eq

	// calc dudx
	nn = CreatedphidxSystemTDSUCS(cUnp, cNXU, cNZU, 0, cUnp, cNXU, cNZU, 0, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDOUCS3(cUnp, cNXU, cNZU, 0, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDUPWIND(cUnp, cNXU, cNZU, 0, cUnp, cNXU, cNZU, 0, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);

	// calc dudz
	nn = CreatedphidzSystemTDSUCS(cUnp, cNXU, cNZU, 0, wnph, (cNXW - 1), cNZW - 1, 2, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDOUCS3(cUnp, cNXU, cNZU, 0, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDUPWIND(cUnp, cNXU, cNZU, 0, wnph, (cNXW - 1), cNZW - 1, 2, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	/*SaveIKCutToCSV("unp", cUnp, cNXU, cNYU, cNZU, cNYU / 2);

	SaveIKCutToCSV("dudx", dphidx, cNXU, cNYU, (cNZU), cNYU / 2);
	SaveIKCutToCSV("dudy", dphidy, cNXU, cNYU, (cNZU), cNYU / 2);
	SaveIKCutToCSV("dudz", dphidz, cNXU, cNYU, (cNZU), cNYU / 2);*/

	// calc d2udx2
	nn = Created2phidx2SystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);


	// calc d2udz2
	nn = Created2phidz2SystemTD(cUnp, cNXU, cNZU, 0, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	// get dudt

	for (int j = 0; j < 1; j++)
	{
#pragma omp parallel for
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(
							cUnp(i, j, k) * dphidx(i, j, k) +
							
							wnph(i, j, k - 1) * dphidz(i, j, k))

						+ (cmu / crho)*(d2phidx2(i, j, k)
							+
							d2phidz2(i, j, k));
				}
			}
		}
	}
}


void SquareCylinder2DStaggeredField::CalcdphidtWmomAdv(Field3D& dphidt,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc)
{
	int nn = 0;

	static Field3DKJI unpt;

	static Field3DIKJ dphidx;
	static Field3DKJI dphidz;
	static Field3DIKJ d2phidx2;
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

	nn = CreateInterpolateZPreIntSystemTD(cUnpc, (cNXU - 1), cNZU, 0, aa, bb, cc, unpt);

	LSSolver::SolveTridiagonalDestructive(unpt, nn, aa, bb, cc);
	
	/*SaveIKCutToCSV("vnpt", vnpt, cNXV, (cNYV - 1), cNZV - 1, cNYP / 2);

	SaveIKCutToCSV("unpt", unpt, cNXU-1, cNYU , cNZU - 1, cNYP / 2);*/

	// w mom eq

	// calc dwdx
	nn = CreatedphidxSystemTDSUCS(cWnp, cNXW, cNZW, 2, unpt, cNXU - 1, cNZU - 1, 0, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDOUCS3(cWnp, cNXW, cNZW, 2, aa, bb, cc, dphidx);
	//nn = CreatedphidxSystemTDUPWIND(cWnp, cNXW, cNZW, 2, unpt, cNXU - 1, cNZU - 1, 0, aa, bb, cc, dphidx);

	LSSolver::SolveTridiagonalDestructive(dphidx, nn, aa, bb, cc);

	// calc dwdz
	nn = CreatedphidzSystemTDSUCS(cWnp, cNXW, cNZW, 2, cWnp, cNXW, cNZW, 2, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDOUCS3(cWnp, cNXW, cNZW, 2, aa, bb, cc, dphidz);
	//nn = CreatedphidzSystemTDUPWIND(cWnp, cNXW, cNZW, 2, cWnp, cNXW, cNZW, 2, aa, bb, cc, dphidz);

	LSSolver::SolveTridiagonalDestructive(dphidz, nn, aa, bb, cc);

	// calc d2udx2
	nn = Created2phidx2SystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, d2phidx2);

	LSSolver::SolveTridiagonalDestructive(d2phidx2, nn, aa, bb, cc);
	
	// calc d2udz2
	nn = Created2phidz2SystemTD(cWnp, cNXW, cNZW, 2, aa, bb, cc, d2phidz2);

	LSSolver::SolveTridiagonalDestructive(d2phidz2, nn, aa, bb, cc);

	//SaveIKCutToCSV("dwdx", dphidx, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("dwdy", dphidy, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("dwdz", dphidz, cNXW, cNYW, cNZW, cNYW / 2);

	//SaveIKCutToCSV("d2wdx2", d2phidx2, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("d2wdy2", d2phidy2, cNXW, cNYW, cNZW, cNYW / 2);
	//SaveIKCutToCSV("d2wdz2", d2phidz2, cNXW, cNYW, cNZW, cNYW / 2);

	// get dwdt
	for (int j = 0; j < 1; j++)
	{

#pragma omp parallel for
		for (int k = 1; k < cNZW - 1; k++)
		{
			for (int i = 1; i < cNXW - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					dphidt(i, j, k) = 0.0;
				}
				else
				{
					dphidt(i, j, k) =
						-1.0*(
							unpt(i - 1, j, k) * dphidx(i, j, k) +
							cWnp(i, j, k) * dphidz(i, j, k))
						+ (cmu / crho)*(d2phidx2(i, j, k)+
							d2phidz2(i, j, k));
				}
			}
		}
	}

	/*SaveIKCutToCSV("dwdt", dphidt, cNXW, cNYW, cNZW, cNYW / 2);*/
}


void SquareCylinder2DStaggeredField::SaveToCSV(const std::string& filename)
{
	std::ofstream arqU(filename + "U.csv");

	for (int k = 0; k < cNZU; k++)
	{
		for (int i = 0; i < cNXU; i++)
		{
			arqU << cUnp(i, 0, k) << "\t";
		}

		arqU << "\n";
	}

	arqU.close();

	std::ofstream arqW(filename + "W.csv");

	for (int k = 0; k < cNZW; k++)
	{
		for (int i = 0; i < cNXW; i++)
		{
			arqW << cWnp(i, 0, k) << "\t";
		}

		arqW << "\n";
	}

	arqW.close();

	std::ofstream arqP(filename + "P.csv");

	for (int k = 0; k < cNZP; k++)
	{
		for (int i = 0; i < cNXP; i++)
		{
			arqP << cPn(i, 0, k) << "\t";
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

int SquareCylinder2DStaggeredField::CreatedphidxUConvSystemTD(
	const Field3D& phiInterpolateX,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = cNXU * 1 * cNZU;

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
		rhs.resize(n2, cNXU, 1, cNZU);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (i == cNXU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if ((i == ciCyInit - 1 || i == ciCyInit - 2 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uipjk = phiInterpolateX(i, j, k);

					const double uimjk = phiInterpolateX(i - 1, j, k);

					rhs(i, j, k) = ((uipjk * uipjk)
						- (uimjk*uimjk)
						) /
						(ch);
				}
				else if (i == 1 || i == 2 || i == cNXU - 2 || i == cNXU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uipjk = phiInterpolateX(i, j, k);

					const double uimjk = phiInterpolateX(i - 1, j, k);

					rhs(i, j, k) = ((uipjk * uipjk)
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

					rhs(i, j, k) = ((uipjk * uipjk)
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

					rhs(i, j, k) = 0.5 / ch * (-3 * uijk * uijk +
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

					rhs(i, j, k) = 0.5 / ch * (3 * uijk*uijk -
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

					rhs(i, j, k) = (1.0 / ch) * ((2 * betaj2 / 3 - 1.0 / 3.0)*uimjk*uimjk -
						(8 * betaj2 / 3 + 0.5)*uijk*uijk
						+ (4 * betaj2 + 1.0)*uipjk*uipjk -
						(8 * betaj2 / 3.0 + 1.0 / 6.0)*uippjk +
						(2 * betaj2 / 3.0*uipppjk*uipppjk));
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

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*uipjk -
						(8 * betajm2 / 3 + 0.5)*uijk
						+ (4 * betajm2 + 1.0)*uimjk -
						(8 * betajm2 / 3.0 + 1.0 / 6.0)*uimmjk +
						(2 * betajm2 / 3.0*uimmmjk));
				}*/
				else
				{

					const double uipjk = phiInterpolateX(i, j, k);
					const double uippjk = phiInterpolateX(i + 1, j, k);

					const double uimjk = phiInterpolateX(i-1, j, k);
					const double uimmjk = phiInterpolateX(i-2, j, k);

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs(i, j, k) = 1.0 / ch * (
						-Fimp / 6.0 * uimmjk * uimmjk
						- Eimp / 1.0*uimjk*uimjk
						+ Eimp / 1.0*uipjk*uipjk
						+ Fimp / 6.0* uippjk * uippjk);*/

					rhs(i, j, k) = bII / (3 * ch)*(uippjk*uippjk - uimmjk)
						+ aII / ch * (uipjk*uipjk - uimjk * uimjk);
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreatedphidzUConvSystemTD(
	const Field3D& phiInterpolateZ,
	const Field3D& W,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = cNXU * 1 * cNZU;

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
		rhs.resize(n2, cNXU, 1, cNZU);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < cNZU; k++)
			for (int i = 0; i < cNXU; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (k == cNZU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if ((k == ckCyInit || k == ckCyInit - 1 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijkp = phiInterpolateZ(i, j, k);
					const double uijkm = phiInterpolateZ(i, j, k-1);

					const double wijkp = W(i, j, k);
					const double wijkm = W(i, j, k-1);

					rhs(i, j, k) = ((uijkp * wijkp)
						- (uijkm * wijkm)) /
						(ch);
				}
				else if (k == 1 || k == 2 || k == cNZU - 2 || k == cNZU - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double uijkp = phiInterpolateZ(i, j, k);
					const double uijkm = phiInterpolateZ(i, j, k-1);

					const double wijkp = W(i, j, k);
					const double wijkm = W(i, j, k-1);

					rhs(i, j, k) = ((uijkp * wijkp)
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

					rhs(i, j, k) = ((uijkp * W[i + ((k) * (cNXW-1)) + ((j)* (cNXW - 1) * cNZW)])
						- (uijkm * W[i + ((k-1) * (cNXW - 1)) + (j  * (cNXW - 1) * cNZW)])) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)] +
						4 * cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
						cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)]);
				}
				else if (k == cNZU - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)] -
						4 * cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] +
						cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * betaj2 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + (k * cNXU) + (j * cNXU * cNZW)]
						+ (4 * betaj2 + 1.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 2) * cNXU) + ((j)* cNXU * cNZW)] +
						(2 * betaj2 / 3.0*cUnp[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k + 3) * cNXU) + ((j)* cNXU * cNZW)]));
				}
				else if (k == cNZU - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cUnp[i + ((k + 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * betajm2 / 3 + 0.5)*cUnp[i + (k * cNXU) + (j * cNXU * cNZU)] * W[i + ((k-1) * cNXU) + (j * cNXU * cNZW)]
						+ (4 * betajm2 + 1.0)*cUnp[i + ((k - 1) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZW)] -
						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cUnp[i + ((k - 2) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZW)] +
						(2 * betajm2 / 3.0*cUnp[i + ((k - 3) * cNXU) + ((j)* cNXU * cNZU)] * W[i + ((k - 4) * cNXU) + ((j)* cNXU * cNZW)]));
				}*/
				else
				{
					const double uijkp = phiInterpolateZ(i, j, k);
					const double uijkpp = phiInterpolateZ(i, j, k+1);

					const double uijkm = phiInterpolateZ(i, j, k-1);
					const double uijkmm = phiInterpolateZ(i, j, k-2);

					const double wijkp = W(i, j, k);
					const double wijkpp = W(i, j, k+1);

					const double wijkm = W(i, j, k-1);
					const double wijkmm = W(i, j, k-2);

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs(i, j, k) = 1.0 / ch * (
						-Fimp / 6.0* uijkmm * wijkmm
						- Eimp / 1.0* uijkm * wijkm
						+ Eimp / 1.0*uijkp * wijkp
						+ Fimp / 6.0*uijkpp * wijkpp);*/

					rhs(i, j, k) = bII / (3 * ch)*(uijkpp*wijkpp - uijkmm * wijkmm)
						+ aII / ch * (uijkp*wijkp - uijkm * wijkm);
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreatedphidxWConvSystemTD(
	const Field3D& phiInterpolateX,
	const Field3D& U,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = cNXW * 1 * cNZW;

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
		rhs.resize(n2, cNXW, 1, cNZW);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (i == cNXW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if ((i == ciCyInit || i == ciCyInit - 1 || i == ciCyEnd - 1 || i == ciCyEnd)
					&& k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wipjk = phiInterpolateX(i, j, k);

					const double wimjk = phiInterpolateX(i-1, j, k);

					const double uipjk = U(i, j, k);
					const double uimjk = U(i-1, j, k);

					rhs(i, j, k) = ((wipjk * uipjk)
						- (wimjk * uimjk)) /
						(ch);
				}
				else if (i == 1 || i == 2 || i == cNXW - 2 || i == cNXW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wipjk = phiInterpolateX(i, j, k);

					const double wimjk = phiInterpolateX(i-1, j, k);

					const double uipjk = U(i, j, k);
					const double uimjk = U(i-1, j, k);

					rhs(i, j, k) = ((wipjk * uipjk)
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

					rhs(i, j, k) = ((wipjk * U[i + ((k)* cNXU) + ((j)* cNXU  * (cNZU-1))])
						- (wimjk * U[i-1 + ((k) * cNXU) + ((j)* cNXU * (cNZU - 1))])) /
						(ch);
				}
				else if (i == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] +
						4 * cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)]);
				}
				else if (i == cNXW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						4 * cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] +
						cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 3 + (k * cNXU) + (j * cNXU * cNZW)]);
				}
				else if (i == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * betaj2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)]
						+ (4 * betaj2 + 1.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 1 + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cWnp[i + 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 2 + (k * cNXU) + (j * cNXU * cNZW)] +
						(2 * betaj2 / 3.0*cWnp[i + 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + 3 + (k * cNXU) + (j * cNXU * cNZW)]));
				}
				else if (i == cNXW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cWnp[i + 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * betajm2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * U[i-1 + (k * cNXU) + (j * cNXU * cNZW)]
						+ (4 * betajm2 + 1.0)*cWnp[i - 1 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 2 + (k * cNXU) + (j * cNXU * cNZW)] -
						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cWnp[i - 2 + (k * cNXW) + (j * cNXW * cNZW)] * U[i -3 + (k * cNXU) + (j * cNXU * cNZW)] +
						(2 * betajm2 / 3.0*cWnp[i - 3 + (k * cNXW) + (j * cNXW * cNZW)] * U[i - 4 + (k * cNXU) + (j * cNXU * cNZW)]));
				}*/
				else
				{
					const double wipjk = phiInterpolateX(i, j, k);
					const double wippjk = phiInterpolateX(i+1, j, k);

					const double wimjk = phiInterpolateX(i-1, j, k);
					const double wimmjk = phiInterpolateX(i-2, j, k);

					const double uipjk = U(i, j, k);
					const double uippjk = U(i+1, j, k);

					const double uimjk = U(i-1, j, k);
					const double uimmjk = U(i-2, j, k);

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*wimmjk * uimmjk -
						Eimp / 1.0*wimjk * uimjk +
						Eimp / 1.0*wipjk * uipjk +
						Fimp / 6.0*wippjk * uippjk);*/

					rhs(i, j, k) = bII / (3 * ch)*(wippjk*uippjk - wimmjk * uimmjk)
						+ aII / ch * (wipjk*uipjk - wimjk * uimjk);
				}
			}

	return n2;
}


int SquareCylinder2DStaggeredField::CreatedphidzWConvSystemTD(
	const Field3D& phiInterpolateZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = cNXW * 1 * cNZW;

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
		rhs.resize(n2, cNXW, 1, cNZW);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < cNZW; k++)
			for (int i = 0; i < cNXW; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if (k == cNZW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;
				}
				else if ((k == ckCyInit - 1 || k == ckCyInit - 2 || k == ckCyEnd - 1 || k == ckCyEnd)
					&& i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wijkp = phiInterpolateZ(i,j ,k);

					const double wijkm = phiInterpolateZ(i, j, k-1);

					rhs(i, j, k) = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				else if (k == 1 || k == 2 || k == cNZW - 2 || k == cNZW - 3)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					const double wijkp = phiInterpolateZ(i, j, k);

					const double wijkm = phiInterpolateZ(i, j, k-1);

					rhs(i, j, k) = ((wijkp*wijkp)
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

					rhs(i, j, k) = ((wijkp*wijkp)
						- (wijkm*wijkm)) /
						(ch);
				}
				else if (k == 0)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (-3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] +
						4 * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == cNZW - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5 / ch * (3 * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] -
						4 * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] +
						cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)]);
				}
				else if (k == 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 1.0 / ch * ((2 * betaj2 / 3 - 1.0 / 3.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * betaj2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNXW) + (j * cNXW * cNZW)]
						+ (4 * betaj2 + 1.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * betaj2 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(2 * betaj2 / 3.0*cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				else if (k == cNZW - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = -1.0 / ch * ((2 * betajm2 / 3 - 1.0 / 3.0)*cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k + 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * betajm2 / 3 + 0.5)*cWnp[i + (k * cNXW) + (j * cNXW * cNZW)] * cWnp[i + (k * cNZW) + (j * cNXW * cNZW)]
						+ (4 * betajm2 + 1.0)*cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW *cNZW)] * cWnp[i + ((k - 1) * cNXW) + ((j)* cNXW * cNZW)] -
						(8 * betajm2 / 3.0 + 1.0 / 6.0)*cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 2) * cNXW) + ((j)* cNXW * cNZW)] +
						(2 * betajm2 / 3.0*cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)] * cWnp[i + ((k - 3) * cNXW) + ((j)* cNXW * cNZW)]));
				}
				*/
				else
				{
					const double wijkp = phiInterpolateZ(i, j, k);
					const double wijkpp = phiInterpolateZ(i, j, k+1);

					const double wijkm = phiInterpolateZ(i, j, k-1);
					const double wijkmm = phiInterpolateZ(i, j, k-2);

					// im
					//aa[index] = Dimp;
					aa[index] = alfaII;

					// i
					bb[index] = 1.0;

					// ip
					//cc[index] = Dimp;
					cc[index] = alfaII;

					/*rhs(i, j, k) = 1.0 / ch * (-Fimp / 6.0*wijkmm * wijkmm -
						Eimp / 1.0*wijkm * wijkm +
						Eimp / 1.0*wijkp * wijkp +
						Fimp / 6.0*wijkpp * wijkpp);*/

					rhs(i, j, k) = bII / (3 * ch)*(wijkpp*wijkpp - wijkmm * wijkmm)
						+ aII / ch * (wijkp*wijkp - wijkm * wijkm);
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreateInterpolateXPreIntSystemTD(
	const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * 1 * NZ;

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
		rhs.resize(n2, NXm, 1, NZ);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				else if (i == NXm - 1 || i == NXm - 2)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				// mantendo condições de contorno ATENÇÃO TEM QUE ESTAR ANTES DO DE BAIXO (ORDEM DOS IFS IMPORTA)
				else if (uvw == 2 && (i == ciCyInit - 1 || i == ciCyEnd - 1) && k >= ckCyInit - 1 && k < ckCyEnd - 1)
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
				else if (uvw == 2 && i >= ciCyInit - 3 && i < ciCyEnd + 2 && k >= ckCyInit - 1 && k < ckCyEnd - 1)
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

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i + 1, j, k)) / 2.0
						+ bI * (phi(i - 1, j, k) + phi(i + 2, j, k)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreateInterpolateZPreIntSystemTD(
	const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * 1 * NZm;

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
		rhs.resize(n2, NX, 1, NZm);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}				
				// MANTENDO COND. CONTORNO ORDEM IMPORTA
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd - 1
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
				else if (uvw == 0 && i >= ciCyInit - 1 && i < ciCyEnd - 1
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

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i, j, k + 1)) / 2.0
						+ bI * (phi(i, j, k - 1) + phi(i, j, k + 2)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreateInterpolateXSystemTD(
	const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * 1 * NZ;

	aa.clear();
	bb.clear();
	cc.clear();
	rhs.cField.clear();

	rhs.cNX = NXm;
	rhs.cNY = 1;
	rhs.cNZ = NZ;

	if (aa.size() != n2)
		aa.resize(n2);

	if (bb.size() != n2)
		bb.resize(n2);

	if (cc.size() != n2)
		cc.resize(n2);

	if (rhs.cField.size() != n2)
		rhs.resize(n2, NXm, 1, NZ);

	for (int j = 0; j < 1; j++)
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
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				else if (i == NXm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw != 0)
						rhs(i, j, k) = 0.0;// 0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				else if (uvw == 0 && (i >= ciCyInit - 3/* || i == ciCyEnd - 1 || i == ciCyEnd ||*/ && i < ciCyEnd + 1)
					&& k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i + 1, j, k));
				}
				else if (uvw == 2 && (i >= ciCyInit - 1 && i <= ciCyEnd - 1)
					&& k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + 1 + (k*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 2 && (i == ciCyInit - 2 || i == ciCyEnd || i == ciCyInit - 3 || i == ciCyEnd + 1)
					&& k >= ckCyInit - 1 && k < ckCyEnd)
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

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i + 1, j, k)) / 2.0
						+ bI * (phi(i - 1, j, k) + phi(i + 2, j, k)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreateInterpolateZSystemTD(
	const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * 1 * NZm;

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
		rhs.resize(n2, NX, 1, NZm);

	for (int j = 0; j < 1; j++)
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
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				else if (k == NZm - 1)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					if (uvw == 1) // só seta igual a zero se for velocidade no y
						rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
					else
						rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				else if (uvw == 0 && (k >= ckCyInit - 1 && k <= ckCyEnd - 1)
					&& i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 0 && (k == ckCyInit - 2 || k == ckCyEnd || k == ckCyInit - 3 || k == ckCyEnd + 1)
					&& i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				else if (uvw == 1 && (k >= ckCyInit - 1 && k == ckCyEnd - 1)
					&& i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.0;//0.5*(phi[i + (k*NX) + (j*NX*NZ)] + phi[i + ((k + 1)*NX) + (j*NX*NZ)]);
				}
				else if (uvw == 1 && (k == ckCyInit - 2 || k == ckCyEnd || k == ckCyInit - 3 || k == ckCyEnd + 1)
					&& i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = 0.5*(phi(i, j, k) + phi(i, j, k + 1));
				}
				else if (uvw == 2 && (k >= ckCyInit - 3 && k < ckCyEnd + 1)
					&& i >= ciCyInit && i < ciCyEnd)
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

					rhs(i, j, k) = aI * (phi(i, j, k) + phi(i, j, k + 1)) / 2.0
						+ bI * (phi(i, j, k - 1) + phi(i, j, k + 2)) / 2.0;
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::CreatedphidxStaggeredSystemTD(const Field3D& phi,
	const int NX, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const
{
	const int NXm = NX - 1;
	const int n2 = NXm * 1 * NZ;

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
		rhs.resize(n2, NXm, 1, NZ);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0 || i == 1 /*|| true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
				}
				else if (i == NXm - 1 || i == NXm - 2)
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

int SquareCylinder2DStaggeredField::CreatedphidzStaggeredSystemTD(const Field3D& phi,
	const int NX, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * 1 * NZm;

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
		rhs.resize(n2, NX, 1, NZm);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 || k == 1/* || true*/)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
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

int SquareCylinder2DStaggeredField::CreatedphidxStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const

{
	const int NXm = NX - 1;
	const int n2 = NXm * 1 * NZ;

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
		rhs.resize(n2, NXm, 1, NZ);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (i == 0 || true)
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

int SquareCylinder2DStaggeredField::CreatedphidzStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NZ,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * 1 * NZm;

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
		rhs.resize(n2, NX, 1, NZm);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (k == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
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

int SquareCylinder2DStaggeredField::CreatedphidxStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs) const

{
	const int NXm = NX - 1;
	const int n2 = NXm * 1 * NZ;

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
		rhs.resize(n2, NXm, 1, NZ);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZ; k++)
			for (int i = 0; i < NXm; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (uvw == 0 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - phi(i, j, k)) / ch;
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


int SquareCylinder2DStaggeredField::CreatedphidzStaggeredSystemTDSimplified(const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs) const
{
	const int NZm = NZ - 1;
	const int n2 = NX * 1* NZm;

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
		rhs.resize(n2, NX, 1, NZm);

	for (int j = 0; j < 1; j++)
		for (int k = 0; k < NZm; k++)
			for (int i = 0; i < NX; i++)
			{
				const int index = rhs.GetIndex(i, j, k);

				if (uvw == 2 || true)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - phi(i, j, k)) / ch;
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

int SquareCylinder2DStaggeredField::Created2phidx2SystemTD(const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DIKJ& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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
						rhs(i, j, k) = (2.0 * -phi(i + 1, j, k) - 5.0*phi(i + 1, j, k) +
							4.0 * phi(i + 2, j, k) - phi(i + 3, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i + 1, j, k) +
							4.0 * phi(i + 2, j, k) - phi(i + 3, j, k)) / (ch*ch);
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
						rhs(i, j, k) = (2.0 * -phi(i - 1, j, k) - 5.0*phi(i - 1, j, k) +
							4.0 * phi(i - 2, j, k) - phi(i - 3, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (2.0 * phi(i, j, k) - 5.0*phi(i - 1, j, k) +
							4.0 * phi(i - 2, j, k) - phi(i - 3, j, k)) / (ch*ch);
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
						rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
							- phi(i, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
							+ phi(i - 1, j, k)) / (ch*ch);
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
							+ phi(i - 1, j, k)) / (ch*ch);
					}
					else
					{
						rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
							+ phi(i - 1, j, k)) / (ch*ch);
					}
				}
				else if (uvw == 0 && i == ciCyInit - 2 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (0.0 - 2.0 *phi(i, j, k)
						+ phi(i - 1, j, k)) / (ch*ch);
				}
				else if (uvw == 0 && i == ciCyInit - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * 0.0 - 5.0*phi(i - 1, j, k) +
						4.0 * phi(i - 2, j, k) - phi(i - 3, j, k)) / (ch*ch);
				}
				else if (uvw == 0 && i == ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
						+ 0.0) / (ch*ch);
				}
				else if (uvw == 0 && i == ciCyEnd - 1 && k >= ckCyInit && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * 0.0 - 5.0*phi(i + 1, j, k) +
						4.0 * phi(i + 2, j, k) - phi(i + 3, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyInit - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i - 1, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyInit && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i - 1, j, k) - 5.0*phi(i - 1, j, k) +
						4.0 * phi(i - 2, j, k) - phi(i - 3, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i + 1, j, k) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				else if (uvw == 2 && i == ciCyEnd - 1 && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i + 1, j, k) - 5.0*phi(i + 1, j, k) +
						4.0 * phi(i + 2, j, k) - phi(i + 3, j, k)) / (ch*ch);
				}
				else
				{
					// im
					aa[index] = alfa2imp;

					// i
					bb[index] = 1.0;

					// ip
					cc[index] = alfa2imp;

					rhs(i, j, k) = (b / (4.0*ch*ch))*(phi(i + 2, j, k) - 2.0*phi(i, j, k) + phi(i - 2, j, k))
						+ (a / (ch*ch))*(phi(i + 1, j, k) - 2.0 * phi(i, j, k) + phi(i - 1, j, k));
				}
			}

	return n2;
}

int SquareCylinder2DStaggeredField::Created2phidz2SystemTD(const Field3D& phi,
	const int NX, const int NZ,
	const int uvw,
	std::vector<double>& aa,
	std::vector<double>& bb,
	std::vector<double>& cc,
	Field3DKJI& rhs)const
{
	const int n2 = NX * 1 * NZ;

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
		rhs.resize(n2, NX, 1, NZ);

	for (int j = 0; j < 1; j++)
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
						rhs(i, j, k) = (2.0 * -phi(i, j, k + 1) - 5.0*phi(i, j, k + 1) +
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
				else if (uvw == 0 && k == ckCyInit - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (-phi(i, j, k) - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyInit && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j, k - 1) - 5.0*phi(i, j, k - 1) +
						4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyEnd - 1 && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * -phi(i, j, k + 1) - 5.0*phi(i, j, k + 1) +
						4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
				}
				else if (uvw == 0 && k == ckCyEnd && i >= ciCyInit - 1 && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						- phi(i, j, k)) / (ch*ch);
				}
				
				else if (uvw == 2 && k == ckCyInit - 2 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (0.0 - 2.0 *phi(i, j, k)
						+ phi(i, j, k - 1)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyInit - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * 0.0 - 5.0*phi(i, j, k - 1) +
						4.0 * phi(i, j, k - 2) - phi(i, j, k - 3)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyEnd - 1 && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (2.0 * 0.0 - 5.0*phi(i, j, k + 1) +
						4.0 * phi(i, j, k + 2) - phi(i, j, k + 3)) / (ch*ch);
				}
				else if (uvw == 2 && k == ckCyEnd && i >= ciCyInit && i < ciCyEnd)
				{
					aa[index] = 0.0;
					cc[index] = 0.0;
					bb[index] = 1.0;

					rhs(i, j, k) = (phi(i, j, k + 1) - 2.0 *phi(i, j, k)
						+ 0.0) / (ch*ch);
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

void SquareCylinder2DStaggeredField::EnforceContinuityVelnp(double dt,
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
		nn = CreatedphidxStaggeredSystemTDSimplified(cUnp, cNXU, cNZU, 0,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);


		nn = CreatedphidzStaggeredSystemTDSimplified(cWnp, cNXW, cNZW, 2,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}
	else
	{
		nn = CreatedphidxStaggeredSystemTD(cUnp, cNXU, cNZU,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);


		nn = CreatedphidzStaggeredSystemTD(cWnp, cNXW, cNZW,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}

	/*SaveIKCutToCSV("dudx", cdphidx, cNXU, (cNYU), cNZU, cNYU / 2);
	SaveIKCutToCSV("dvdy", cdphidy, cNXV, (cNYV), cNZV, cNYV / 2);
	SaveIKCutToCSV("dwdz", cdphidz, cNXW, (cNYW), cNZW, cNYW / 2);

	SaveIKCutToCSV("Unp", cUnp, cNXU, (cNYU), cNZU, cNYU / 2);
	SaveIKCutToCSV("Vnp", cVnp, cNXV, (cNYV), cNZV, cNYV / 2);
	SaveIKCutToCSV("Wnp", cWnp, cNXW, (cNYW), cNZW, cNYW / 2);*/

	nn = CreateSparsePressureEquation(cdphidx, cdphidz,
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

	if (show)
		std::cout << "Pressure eq solved its " << nitPEqt << " err " << errOutPEqt << std::endl;

	//SaveIKCutToCSV("press", cPn, cNXP, (cNYP), cNZP, cNYP / 2);

	// project p into uvw

	if (true)
	{
		nn = CreatedphidxStaggeredSystemTDSimplified(cPn, cNXP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);
		
		nn = CreatedphidzStaggeredSystemTDSimplified(cPn, cNXP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}
	else
	{
		// set dp/dn = 0
		for (int j = 0; j < 1; j++)
			for (int k = 0; k < cNZP; k++)
			{
				cPn(0, j, k) = cPn(1, j, k);
				cPn(cNXP - 1, j, k) = cPn(cNXP - 2, j, k);
			}

		nn = CreatedphidxStaggeredSystemTD(cPn, cNXP, cNZP,
			aa, bb, cc, cdphidx);

		lssolver.SolveTridiagonalDestructive(cdphidx, nn, aa, bb, cc);
		
		// set dp/dn = 0
		for (int j = 0; j < 1; j++)
			for (int i = 0; i < cNXP; i++)
			{
				cPn(i, j, 0) = cPn(i, j, 1);
				cPn(i, j, cNZP - 1) = cPn(i, j, cNZP - 2);
			}

		nn = CreatedphidzStaggeredSystemTD(cPn, cNXP, cNZP,
			aa, bb, cc, cdphidz);

		lssolver.SolveTridiagonalDestructive(cdphidz, nn, aa, bb, cc);
	}

	//SaveIKCutToCSV("dPdx", cdphidx, cNXP, (cNYP), cNZP, cNYP / 2);
	//SaveIKCutToCSV("dPdy", cdphidy, cNXP, (cNYP), cNZP, cNYP / 2);
	//SaveIKCutToCSV("dPdz", cdphidz, cNXP, (cNYP), cNZP, cNYP / 2);


	for (int j = 0; j < 1; j++)
	{
		for (int k = 1; k < cNZU - 1; k++)
		{
			for (int i = 1; i < cNXU - 1; i++)
			{
				if (i >= ciCyInit - 1 && i < ciCyEnd && k >= ckCyInit && k < ckCyEnd)
				{
					continue;
				}
				cUnp(i, j, k) -=
					(dt / crho)*(cdphidx(i, j, k));

			}
		}
	}
	

	for (int j = 0; j < 1; j++)
	{
		for (int k = 1; k < cNZW - 1; k++)
		{
			for (int i = 1; i < cNXW - 1; i++)
			{
				if (i >= ciCyInit && i < ciCyEnd && k >= ckCyInit - 1 && k < ckCyEnd)
				{
					continue;
				}
				cWnp(i, j, k) -=
					(dt / crho)*(cdphidz(i, j, k));

			}
		}
	}
}

void SquareCylinder2DStaggeredField::SaveField(const std::string& baseFileName) const
{
	const int64_t NY = 1;
	std::ofstream uFile(baseFileName + "U.dat", std::ios::binary);

	uFile.write(reinterpret_cast<const char*>(&cNXU), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&NY), sizeof(int64_t));
	uFile.write(reinterpret_cast<const char*>(&cNZU), sizeof(int64_t));

	uFile.write(reinterpret_cast<const char*>(cUnp.data()), cUnp.size() * sizeof(double));

	uFile.close();

	std::ofstream wFile(baseFileName + "W.dat", std::ios::binary);

	wFile.write(reinterpret_cast<const char*>(&cNXW), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&NY), sizeof(int64_t));
	wFile.write(reinterpret_cast<const char*>(&cNZW), sizeof(int64_t));

	wFile.write(reinterpret_cast<const char*>(cWnp.data()), cWnp.size() * sizeof(double));

	wFile.close();

	std::ofstream pFile(baseFileName + "P.dat", std::ios::binary);

	pFile.write(reinterpret_cast<const char*>(&cNXP), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&NY), sizeof(int64_t));
	pFile.write(reinterpret_cast<const char*>(&cNZP), sizeof(int64_t));

	pFile.write(reinterpret_cast<const char*>(cPn.data()), cPn.size() * sizeof(double));

	pFile.close();
}

void SquareCylinder2DStaggeredField::ReadField(const std::string& basefilename)
{
	int64_t NY = 1;
	std::ifstream uFile(basefilename + "U.dat", std::ios::binary);

	uFile.read(reinterpret_cast<char*>(&cNXU), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&NY), sizeof(int64_t));
	uFile.read(reinterpret_cast<char*>(&cNZU), sizeof(int64_t));

	cUn.clear();
	cUn.resize(cNXU * NY * cNZU, cNXU, NY, cNZU);

	uFile.read(reinterpret_cast<char*>(cUn.data()), cUn.size() * sizeof(double));

	uFile.close();

	cUnp = cUn;

	std::ifstream wFile(basefilename + "W.dat", std::ios::binary);

	wFile.read(reinterpret_cast<char*>(&cNXW), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&NY), sizeof(int64_t));
	wFile.read(reinterpret_cast<char*>(&cNZW), sizeof(int64_t));

	cWn.clear();
	cWn.resize(cNXW * NY * cNZW, cNXW, NY, cNZW);

	wFile.read(reinterpret_cast<char*>(cWn.data()), cWn.size() * sizeof(double));

	wFile.close();

	cWnp = cWn;

	std::ifstream pFile(basefilename + "P.dat", std::ios::binary);

	pFile.read(reinterpret_cast<char*>(&cNXP), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&NY), sizeof(int64_t));
	pFile.read(reinterpret_cast<char*>(&cNZP), sizeof(int64_t));

	cPn.clear();
	cPn.resize(cNXP * NY * cNZP, cNXP, NY, cNZP);

	pFile.read(reinterpret_cast<char*>(cPn.data()), cPn.size() * sizeof(double));

	pFile.close();
}

void SquareCylinder2DStaggeredField::InterpolateUWInRespDirections(
	std::vector<double>& aa, std::vector<double>& bb, std::vector<double>& cc,
	std::vector<double>& aa1, std::vector<double>& bb1, std::vector<double>& cc1,
	std::vector<double>& aa2, std::vector<double>& bb2, std::vector<double>& cc2)
{
	SquareCylinder2DStaggeredField& classe = *this;


	std::thread thread1([&aa, &bb, &cc, &classe]()
		{
			int nn = 0;
			// interpolate u in x
			nn = classe.CreateInterpolateXSystemTD(classe.cUnp, classe.cNXU, classe.cNZU, 0,
				aa, bb, cc, classe.cUnpc);

			LSSolver::SolveTridiagonalDestructive(classe.cUnpc, nn, aa, bb, cc);
		});
	
	int nn = 0;
	// interpolate w in z
	nn = CreateInterpolateZSystemTD(cWnp, cNXW, cNZW, 2,
		aa2, bb2, cc2, cWnpc);

	LSSolver::SolveTridiagonalDestructive(cWnpc, nn, aa2, bb2, cc2);

	thread1.join();

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