#pragma once
#include <vector>
#include <fstream>

class FDColocated3DField
{
public:
	typedef std::vector<double> field3D;

	const double Dimp = 0.3793894912;
	const double Eimp = 1.57557379;
	const double Fimp = 0.183205192;

	const double alfa2imp = 0.26;
	const double a = (4.0 / 3.0)*(1.0 - alfa2imp);
	const double b = (1.0 / 3.0)*(10.0*alfa2imp - 1.0);

	const double alfaI = 3.0 / 10.0; //0.42;
	const double aI = (9.0 + 10.0*alfaI) / 8.0;
	const double bI = (6 * alfaI - 1.0) / 8.0;

	const double alfaII = 9.0 / 62.0;
	const double aII = (9.0 - 6.0*alfaII) / 8.0;
	const double bII = (22.0*alfaII - 1.0) / 8.0;

public:
	FDColocated3DField(double rho, double mu, double h, double L, double H, double W, double freestream);
	FDColocated3DField() = delete;
	~FDColocated3DField();

	void SetBCFlatPlate();
	void UpdateBCFlatPlate();

	void SaveIKCutToCSV(const std::string& filename, int j);

	void SaveJKCutToCSV(const std::string& filename, int i);

	int CreatedphidxSystemTDUNUSED(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa, 
		std::vector<double>& bb, 
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidySystemTDUNUSED(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidzSystemTDUNUSED(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa, 
		std::vector<double>& bb, 
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidxStaggeredSystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidyStaggeredSystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidzStaggeredSystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidxStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidyStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidzStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int Created2phidx2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int Created2phidy2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int Created2phidz2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int Created2phidx2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int Created2phidy2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int Created2phidz2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidxUConvSystemTD(
		std::vector<double>& aa, 
		std::vector<double>& bb, 
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidyUConvSystemTD(
		const field3D& V,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidzUConvSystemTD(
		const field3D& W,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidxVConvSystemTD(
		const field3D& U,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidyVConvSystemTD(
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidzVConvSystemTD(
		const field3D& W,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidxWConvSystemTD(
		const field3D& U,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidyWConvSystemTD(
		const field3D& V,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreatedphidzWConvSystemTD(
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreateInterpolateXSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreateInterpolateYSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreateInterpolateZSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<double>& rhs) const;

	int CreateSparsePressureEquation(const field3D& dUdx,
		const field3D& dVdy,
		const field3D& dWdz,
		const double dt,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs) const;

	void CalcdphidtUmom(field3D& dphidt, 
		class LSSolver& lssolver, 
		std::vector<double>& aa, 
		std::vector<double>& bb, 
		std::vector<double>& cc);

	void CalcdphidtVmom(field3D& dphidt,
		LSSolver& lssolver,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc);

	void CalcdphidtWmom(field3D& dphidt,
		LSSolver& lssolver,
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc);

	void EnforceContinuityVelnp(double dt, 
		std::vector<double>& aa,
		std::vector<double>& bb,
		std::vector<double>& cc,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<double>& val,
		std::vector<double>& rhs,
		LSSolver& lssolver);
	
public:
	const double crho;
	const double cmu;
	const double ch;

	const double cL;
	const double cH;
	const double cW;

	const int cNXU;
	const int cNYU;
	const int cNZU;

	const int cNXV;
	const int cNYV;
	const int cNZV;

	const int cNXW;
	const int cNYW;
	const int cNZW;

	const int cNXP;
	const int cNYP;
	const int cNZP;

	const double cUfreestream;

	bool cFirstIteration;

	field3D cUn;
	field3D cVn;
	field3D cWn;
	field3D cPn;

	field3D cUnp;
	field3D cVnp;
	field3D cWnp;
	field3D cPnp;

	field3D cPhi;

	field3D cdphidx;
	field3D cdphidy;
	field3D cdphidz;

	field3D cphiInterpolateX;
	field3D cphiInterpolateY;
	field3D cphiInterpolateZ;

	field3D cd2phidx2;
	field3D cd2phidy2;
	field3D cd2phidz2;

	field3D cdUdtn;
	field3D cdUdt1;
	field3D cdUdt2;
	field3D cdUdt3;

	field3D cdVdtn;
	field3D cdVdt1;
	field3D cdVdt2;
	field3D cdVdt3;

	field3D cdWdtn;
	field3D cdWdt1;
	field3D cdWdt2;
	field3D cdWdt3;
};