#pragma once
#include <vector>
#include <fstream>

class SquareCylinder3DStaggeredFieldFFast
{
public:
	typedef std::vector<float> field3D;

	const float Dimp = 0.3793894912;
	const float Eimp = 1.57557379;
	const float Fimp = 0.183205192;

	const float alfa2imp = 0.26;
	const float a = (4.0 / 3.0)*(1.0 - alfa2imp);
	const float b = (1.0 / 3.0)*(10.0*alfa2imp - 1.0);

	const float alfaI = 0.42;//3.0 / 10.0;
	const float aI = (9.0 + 10.0*alfaI) / 8.0;
	const float bI = (6 * alfaI - 1.0) / 8.0;

	const float alfaII = 0.216;//9.0 / 62.0;
	const float aII = (9.0 - 6.0*alfaII) / 8.0;
	const float bII = (22.0*alfaII - 1.0) / 8.0;

public:
	SquareCylinder3DStaggeredFieldFFast(float rho, float mu, float h, float l, float L, float H, float W, float freestream);
	SquareCylinder3DStaggeredFieldFFast() = delete;
	~SquareCylinder3DStaggeredFieldFFast();

	void SetBCFlatPlate();
	void UpdateBCFlatPlate();

	void SaveIKCutToCSV(const std::string& filename, int j);

	void SaveJKCutToCSV(const std::string& filename, int i);

	void SaveIJCutToCSV(const std::string& filename, int k);

	int CreatedphidxSystemTDUNUSED(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidySystemTDUNUSED(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzSystemTDUNUSED(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxStaggeredSystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyStaggeredSystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzStaggeredSystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzStaggeredSystemTDSimplified(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int Created2phidx2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int Created2phidy2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int Created2phidz2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int Created2phidx2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int Created2phidy2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int Created2phidz2SystemTD(const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxUConvSystemTD(
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyUConvSystemTD(
		const field3D& V,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzUConvSystemTD(
		const field3D& W,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxVConvSystemTD(
		const field3D& U,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyVConvSystemTD(
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzVConvSystemTD(
		const field3D& W,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxWConvSystemTD(
		const field3D& U,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyWConvSystemTD(
		const field3D& V,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzWConvSystemTD(
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateInterpolateXSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateInterpolateYSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateInterpolateZSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateInterpolateXSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateInterpolateYSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateInterpolateZSystemTD(
		const field3D& phi,
		const int NX, const int NY, const int NZ,
		const int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreateSparsePressureEquation(const field3D& dUdx,
		const field3D& dVdy,
		const field3D& dWdz,
		const float dt,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<float>& val,
		std::vector<float>& rhs) const;

	void CalcdphidtUmom(field3D& dphidt,
		class LSSolverF& lssolver,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtVmom(field3D& dphidt,
		LSSolverF& lssolver,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtWmom(field3D& dphidt,
		LSSolverF& lssolver,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void EnforceContinuityVelnp(float dt,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<int>& ptr,
		std::vector<int>& col,
		std::vector<float>& val,
		std::vector<float>& rhs,
		LSSolverF& lssolver);

	void SaveField(const std::string& baseFileName) const;

	void ReadField(const std::string& basefilename);

public:
	const float crho;
	const float cmu;
	const float ch;

	const float cL;
	const float cH;
	const float cW;

	const float cl;

	const int cDistRelFrenteRect;

	int cNXU;
	int cNYU;
	int cNZU;

	int cNXV;
	int cNYV;
	int cNZV;

	int cNXW;
	int cNYW;
	int cNZW;

	int cNXP;
	int cNYP;
	int cNZP;

	int cNXCy;
	int cNYCy;
	int cNZCy;

	int ciCyInit;
	int cjCyInit;
	int ckCyInit;

	int ciCyEnd;
	int cjCyEnd;
	int ckCyEnd;

	float cUfreestream;

	bool cFirstIteration;

	field3D cUn;
	field3D cVn;
	field3D cWn;
	field3D cPn;

	field3D cUnp;
	field3D cVnp;
	field3D cWnp;
	/*field3D cPnp;*/

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