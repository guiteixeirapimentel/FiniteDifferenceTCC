#pragma once
#include <iostream>

#include "FDColocated3DField.h"
#include "LSSolver.h"

#include "Timer.h"

void LidDriven()
{
	constexpr double h = 0.01;
	
	constexpr double L = 1.0;
	constexpr double H = 1.0;
	constexpr double W = 1.0;

	constexpr int NX = int(L / h + 2);
	constexpr int NY = int(H / h + 2);
	constexpr int NZ = int(W / h + 2);

	constexpr int NN = int(NX * NY * NZ);

	constexpr double ufree = 1.0;
	constexpr double vfree = 0.0;
	constexpr double wfree = 0.0;

	constexpr double dt = h / 20;

	constexpr double cfl = ufree * dt / h;

	constexpr double rho = 1.0;
	constexpr double mu = 1e-2;

	constexpr double tfinal = 15.0;

	constexpr int NT = int(tfinal / dt) + 1;

	constexpr double CellRe = rho * h*ufree / mu;

	constexpr double Re = rho * L*ufree / mu;

	constexpr int saveEvery = 63;

	constexpr double dtSave = saveEvery * dt;

	constexpr double conditionForStability = 4 * dt / (Re*h*h);

	std::vector<double> aa, bb, cc;
	std::vector<double> rhs;
	std::vector<double> val;
	std::vector<int> ptr, col;

	LSSolver lssolver;

	FDColocated3DField fdcolfield(rho, mu, h, L, H, W, ufree);

	fdcolfield.SetBCFlatPlate();

	//fdcolfield.UpdateBCFlatPlate();

	fdcolfield.SaveIKCutToCSV("BC", (h / h + 1) / 2);

	double t = 0.0;

	Timer timer;

	const int idttimestep = timer.SetTimer("Complete timestep");

	const int idtuvwmom1 = timer.SetTimer("uvwmom 1 calculation");
	const int idtpress1 = timer.SetTimer("pressure 1 calculation and continuity");

	const int idtuvwmom2 = timer.SetTimer("uvwmom 2 calculation");
	const int idtpress2 = timer.SetTimer("pressure 2 calculation and continuity");

	const int idtuvwmom3 = timer.SetTimer("uvwmom 3 calculation");
	const int idtpress3 = timer.SetTimer("pressure 3 calculation and continuity");

	const int idtuvwmom4 = timer.SetTimer("uvwmom 4 calculation");
	const int idtpress4 = timer.SetTimer("pressure 4 calculation and continuity");

	for (int it = 0; it < NT; it++)
	{
		int nn = 0;

		//fdcolfield.UpdateBCFlatPlate();

		timer.Tick(idttimestep);

		// momentum eqts.

		timer.Tick(idtuvwmom1);

		fdcolfield.CalcdphidtUmom(fdcolfield.cdUdtn, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtVmom(fdcolfield.cdVdtn, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtWmom(fdcolfield.cdWdtn, lssolver, aa, bb, cc);


		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/

					fdcolfield.cUnp[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] =
						fdcolfield.cUn[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
						(dt / 2.0)*fdcolfield.cdUdtn[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)];
				}
			}
		}

		for (int j = 1; j < fdcolfield.cNYV - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZV - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXV - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cVnp[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] =
						fdcolfield.cVn[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
						(dt / 2.0)*fdcolfield.cdVdtn[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)];
				}
			}
		}

		for (int j = 1; j < fdcolfield.cNYW - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cWnp[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] =
						fdcolfield.cWn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
						(dt / 2.0)*fdcolfield.cdWdtn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)];
				}
			}
		}


		timer.Tock(idtuvwmom1);

		// calculate pressure and enforce continuity


		timer.Tick(idtpress1);

		fdcolfield.EnforceContinuityVelnp(dt/2.0, aa, bb, cc, ptr, col, val, rhs, lssolver);
		
		timer.Tock(idtpress1);

		timer.Tick(idtuvwmom2);

		fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt1, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtVmom(fdcolfield.cdVdt1, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt1, lssolver, aa, bb, cc);

		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cUnp[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] =
						fdcolfield.cUn[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
						(dt / 2.0)*fdcolfield.cdUdt1[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)];
				}
			}
		}

		for (int j = 1; j < fdcolfield.cNYV - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZV - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXV - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cVnp[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] =
						fdcolfield.cVn[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
						(dt / 2.0)*fdcolfield.cdVdt1[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)];
				}
			}
		}

		for (int j = 1; j < fdcolfield.cNYW - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cWnp[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] =
						fdcolfield.cWn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
						(dt / 2.0)*fdcolfield.cdWdt1[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)];
				}
			}
		}

		timer.Tock(idtuvwmom2);


		// calculate pressure and enforce continuity

		timer.Tick(idtpress2);

		fdcolfield.EnforceContinuityVelnp(dt/2.0, aa, bb, cc, ptr, col, val, rhs, lssolver);

		timer.Tock(idtpress2);

		timer.Tick(idtuvwmom3);

		fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt2, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtVmom(fdcolfield.cdVdt2, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt2, lssolver, aa, bb, cc);

		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cUnp[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] =
						fdcolfield.cUn[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
						(dt)*fdcolfield.cdUdt2[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)];
				}
			}
		}


		for (int j = 1; j < fdcolfield.cNYV - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZV - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXV - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cVnp[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] =
						fdcolfield.cVn[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
						(dt)*fdcolfield.cdVdt2[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)];
				}
			}
		}

		for (int j = 1; j < fdcolfield.cNYW - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cWnp[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] =
						fdcolfield.cWn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
						(dt)*fdcolfield.cdWdt2[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)];
				}
			}
		}

		timer.Tock(idtuvwmom3);

		// calculate pressure and enforce continuity

		timer.Tick(idtpress3);

		fdcolfield.EnforceContinuityVelnp(dt, aa, bb, cc, ptr, col, val, rhs, lssolver);

		timer.Tock(idtpress3);

		timer.Tick(idtuvwmom4);

		fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt3, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtVmom(fdcolfield.cdVdt3, lssolver, aa, bb, cc);

		fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt3, lssolver, aa, bb, cc);


		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cUnp[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] =
						fdcolfield.cUn[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
						(dt / 6.0)*(
							fdcolfield.cdUdtn[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
							2.0 * fdcolfield.cdUdt1[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
							2.0 * fdcolfield.cdUdt2[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)] +
							fdcolfield.cdUdt3[i + (k*fdcolfield.cNXU) + (j*fdcolfield.cNXU*fdcolfield.cNZU)]
							);
				}
			}
		}

		for (int j = 1; j < fdcolfield.cNYV - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZV - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXV - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cVnp[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] =
						fdcolfield.cVn[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
						(dt / 6.0)*(
							fdcolfield.cdVdtn[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
							2.0 * fdcolfield.cdVdt1[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
							2.0 * fdcolfield.cdVdt2[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)] +
							fdcolfield.cdVdt3[i + (k*fdcolfield.cNXV) + (j*fdcolfield.cNXV*fdcolfield.cNZV)]
							);
				}
			}
		}

		timer.Tock(idtuvwmom4);

		for (int j = 1; j < fdcolfield.cNYW - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					/*if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}*/
					fdcolfield.cWnp[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] =
						fdcolfield.cWn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
						(dt / 6.0)*(
							fdcolfield.cdWdtn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
							2.0 * fdcolfield.cdWdt1[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
							2.0 * fdcolfield.cdWdt2[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
							fdcolfield.cdWdt3[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)]
							);
				}
			}
		}


		// calculate pressure and enforce continuity

		timer.Tick(idtpress4);

		fdcolfield.EnforceContinuityVelnp(dt, aa, bb, cc, ptr, col, val, rhs, lssolver);

		timer.Tock(idtpress4);

		// START SAVE TO FILE

		{
			nn = fdcolfield.CreatedphidxStaggeredSystemTDSimplified(fdcolfield.cUnp, fdcolfield.cNXU, fdcolfield.cNYU, fdcolfield.cNZU,
				aa, bb, cc, fdcolfield.cdphidx);

			lssolver.SolveTridiagonalDestructive(fdcolfield.cdphidx, nn, aa, bb, cc);

			nn = fdcolfield.CreatedphidyStaggeredSystemTDSimplified(fdcolfield.cVnp, fdcolfield.cNXV, fdcolfield.cNYV, fdcolfield.cNZV,
				aa, bb, cc, fdcolfield.cdphidy);

			lssolver.SolveTridiagonalDestructive(fdcolfield.cdphidy, nn, aa, bb, cc);

			nn = fdcolfield.CreatedphidzStaggeredSystemTDSimplified(fdcolfield.cWnp, fdcolfield.cNXW, fdcolfield.cNYW, fdcolfield.cNZW,
				aa, bb, cc, fdcolfield.cdphidz);

			lssolver.SolveTridiagonalDestructive(fdcolfield.cdphidz, nn, aa, bb, cc);
		}

		if ((it % saveEvery) == 0)
		{
			fdcolfield.SaveIKCutToCSV("flatplateturb/save" + std::to_string(t), int(W / h) / 2);

			fdcolfield.SaveJKCutToCSV("flatplateturb/save" + std::to_string(t), int(L / (h) / 2));
		}

		double maxDiv = 0.0;

		int maxDivi, maxDivj, maxDivk;
		maxDivi = maxDivj = maxDivk = 0;

		for (int j = 1; j < fdcolfield.cNYP - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZP - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXP - 2; i++)
				{
					const double div = fdcolfield.cdphidx[i - 1 + ((k)*(fdcolfield.cNXU - 1)) + ((j) * (fdcolfield.cNXU - 1) * fdcolfield.cNZU)] +
						fdcolfield.cdphidy[i + ((k)*fdcolfield.cNXV) + ((j - 1) * fdcolfield.cNXV * fdcolfield.cNZV)] +
						fdcolfield.cdphidz[i + ((k - 1)*fdcolfield.cNXW) + ((j)* fdcolfield.cNXW * (fdcolfield.cNZW - 1))];

					if (std::fabs(div) > maxDiv)
					{
						maxDiv = div;

						maxDivi = i;
						maxDivj = j;
						maxDivk = k;

					}
				}
			}
		}

		std::cout << "max div" << maxDiv << " i " << maxDivi << " j " << maxDivj << " k " << maxDivk << std::endl;
		std::cout << "tempo " << t << std::endl;

		timer.Tock(idttimestep);

		timer.WriteToCoutAllTickTock();

		// END SAVE TO FILE

		// set old to new
		fdcolfield.cUn = fdcolfield.cUnp;
		fdcolfield.cVn = fdcolfield.cVnp;
		fdcolfield.cWn = fdcolfield.cWnp;

		t += dt;
	}
}