#pragma once
#include <vector>
#include <fstream>

class SquareCylinder3DStaggeredFieldF
{
public:
	typedef std::vector<float> field3D;

	const float Dimp = 0.3793894912;
	const float Eimp = 1.57557379;
	const float Fimp = 0.183205192;
	const float etastar = -2.0;

	const float qm2 = (-Fimp / 4.0) + (etastar / 300.0);
	const float qm1 = (-Eimp / 2.0) + (etastar / 30.0);
	const float qp1 = (Eimp / 2.0) + (etastar / 30.0);
	const float qp2 = (Fimp / 4.0) + (etastar / 300.0);
	const float q0 = (-11.0*etastar) / 150.0;
	const float pip = Dimp + (etastar / 60.0);
	const float pim = Dimp - (etastar / 60.0);

	const float betaj2 = -0.025; //-0.025;
	const float betajm2 = 0.09;  //0.09;

	const float alfa2imp = 2.0 / 11.0;//0.26;
	const float a = (4.0 / 3.0)*(1.0 - alfa2imp);
	const float b = (1.0 / 3.0)*(10.0*alfa2imp - 1.0);

	const float alfaI = 3.0 / 10.0; // 0.42;
	const float aI = (9.0 + 10.0*alfaI) / 8.0;
	const float bI = (6 * alfaI - 1.0) / 8.0;

	const float alfaII = 0.216; //9.0 / 62.0;
	const float aII = (9.0 - 6.0*alfaII) / 8.0;
	const float bII = (22.0*alfaII - 1.0) / 8.0;

	const float velPosSUCSalfam1 = 1.0/2.0;
	const float velPosSUCSalfap1 = 1.0 / 6.0;
	const float velPosSUCScm2 = -1.0 / 18.0;
	const float velPosSUCScm1 = -1.0;
	const float velPosSUCSc0 = 1.0 / 2.0;
	const float velPosSUCSc1 = 5.0 / 9.0;

	const float velNegSUCSalfam1 = 1.0 / 6.0;
	const float velNegSUCSalfap1 = 1.0 / 2.0;
	const float velNegSUCSc0 = -1.0 / 2.0;
	const float velNegSUCSc1 = 1.0;
	const float velNegSUCSc2 = 1.0 / 18.0;

public:
	SquareCylinder3DStaggeredFieldF(float rho, float mu, float h, float l, float L, float H, float W, float freestream);
	SquareCylinder3DStaggeredFieldF() = delete;
	~SquareCylinder3DStaggeredFieldF();

	void SetBCFlatPlate();
	void UpdateBCFlatPlate();

	void SaveIKCutToCSV(const std::string& filename, int j);

	void SaveJKCutToCSV(const std::string& filename, int i);

	void SaveIJCutToCSV(const std::string& filename, int k);

	int CreatedphidxSystemTDOUCS3(const field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidySystemTDOUCS3(const field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzSystemTDOUCS3(const field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxSystemTDSUCS(const field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidySystemTDSUCS(const field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzSystemTDSUCS(const field3D& phi,
		const int NX, const int NY, const int NZ,
		int uvw,
		const field3D& convVel,
		const int NXc, const int NYc, const int NZc,
		int uvwc,
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
		const std::vector<float>& phiInterpolateX,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyUConvSystemTD(
		const std::vector<float>& phiInterpolateY,
		const field3D& V,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzUConvSystemTD(
		const std::vector<float>& phiInterpolateZ,
		const field3D& W,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxVConvSystemTD(
		const std::vector<float>& phiInterpolateX,
		const field3D& U,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyVConvSystemTD(
		const std::vector<float>& phiInterpolateY,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzVConvSystemTD(
		const std::vector<float>& phiInterpolateZ,
		const field3D& W,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidxWConvSystemTD(
		const std::vector<float>& phiInterpolateX,
		const field3D& U,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidyWConvSystemTD(
		const std::vector<float>& phiInterpolateY,
		const field3D& V,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc,
		std::vector<float>& rhs) const;

	int CreatedphidzWConvSystemTD(
		const std::vector<float>& phiInterpolateZ,
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

	void CalcdphidtUmomDiv(field3D& dphidt,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtVmomDiv(field3D& dphidt,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtWmomDiv(field3D& dphidt,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtUmomAdv(field3D& dphidt,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtVmomAdv(field3D& dphidt,
		std::vector<float>& aa,
		std::vector<float>& bb,
		std::vector<float>& cc);

	void CalcdphidtWmomAdv(field3D& dphidt,
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
		class LSSolverF& lssolver);

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

	field3D cdPUdt;
	field3D cdPVdt;
	field3D cdPWdt;
	
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