#pragma once
#include <iostream>
#include <thread>

#include "SquareCylinder2DStaggeredField.h"
#include "LSSolver.h"

#include "Timer.h"

void SquareRect2D()
{
	constexpr double h = 0.0020;//02;

	constexpr double l = 0.02;

	constexpr double L = 26.0*l;
	constexpr double H = 15.0*l;

	constexpr int NX = int(L / h + 2);
	constexpr int NZ = int(H / h + 2);

	constexpr int NN = int(NX * 1 * NZ);

	constexpr double ufree = 0.15;
	constexpr double vfree = 0.0;
	constexpr double wfree = 0.0;

	constexpr double dt = h / 2;

	constexpr double cfl = 4 * ufree * dt / h;

	constexpr double rho = 1.0;
	constexpr double mu = 3.0e-5;

	constexpr double tfinal = 300.0;

	constexpr int NT = int(tfinal / dt) + 1;

	constexpr double CellRe = rho * h*ufree / mu;

	constexpr double Re = rho * l*ufree / mu;

	constexpr int saveEvery = 366 * (0.001 / dt);

	constexpr double dtSave = saveEvery * dt;

	constexpr double conditionForStability = 4 * dt / (Re*h*h);

	std::vector<double> aa, bb, cc;
	Field3DIKJ rhs;
	std::vector<double> val;
	std::vector<int> ptr, col;

	LSSolver lssolver;

	SquareCylinder2DStaggeredField fdcolfield(rho, mu, h, l, L, H, ufree);

	fdcolfield.SetBCFlatPlate();

	fdcolfield.SaveToCSV("BC");

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

		// Set inlet pertubation
		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				fdcolfield.cUnp(1, j, k) =
					fdcolfield.cUfreestream + (sin(t + (k * fdcolfield.ch)) * fdcolfield.cUfreestream * 0.01);

			}
		}

		fdcolfield.UpdateBCFlatPlate();

		timer.Tick(idttimestep);

		// momentum eqts.

		timer.Tick(idtuvwmom1);

		fdcolfield.CalcdphidtUmomDiv(fdcolfield.cdUdtn, aa, bb, cc);
		
		fdcolfield.CalcdphidtWmomDiv(fdcolfield.cdWdtn, aa, bb, cc);


		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd)
					{
						continue;
					}

					fdcolfield.cUnp(i, j, k) =
						fdcolfield.cUn(i, j, k) +
						(dt / 2.0)*fdcolfield.cdUdtn(i, j, k);
				}
			}
		}
		
		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cWnp(i, j, k) =
						fdcolfield.cWn(i, j, k) +
						(dt / 2.0)*fdcolfield.cdWdtn(i, j, k);
				}
			}
		}


		timer.Tock(idtuvwmom1);

		// calculate pressure and enforce continuity

		//fdcolfield.SaveIKCutToCSV("flatplateturb/prepress" + std::to_string(t), int(W / h) / 2);

		timer.Tick(idtpress1);

		//fdcolfield.EnforceContinuityVelnp(dt / 2.0, aa, bb, cc, ptr, col, val, rhs, true, lssolver);

		timer.Tock(idtpress1);

		timer.Tick(idtuvwmom2);

		fdcolfield.CalcdphidtUmomDiv(fdcolfield.cdUdt1, aa, bb, cc);
		
		fdcolfield.CalcdphidtWmomDiv(fdcolfield.cdWdt1, aa, bb, cc);

		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cUnp(i, j, k) =
						fdcolfield.cUn(i, j, k) +
						(dt / 2.0)*fdcolfield.cdUdt1(i, j, k);
				}
			}
		}

		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cWnp(i, j, k) =
						fdcolfield.cWn(i, j, k) +
						(dt / 2.0)*fdcolfield.cdWdt1(i, j, k);
				}
			}
		}

		timer.Tock(idtuvwmom2);


		//fdcolfield.SaveJKCutToCSV("flatplateturb/prepress" + std::to_string(t), int(L / h) / 2);

		// calculate pressure and enforce continuity

		timer.Tick(idtpress2);


		fdcolfield.EnforceContinuityVelnp(dt / 2.0, aa, bb, cc, ptr, col, val, rhs, true, lssolver);

		//fdcolfield.SaveJKCutToCSV("flatplateturb/premom" + std::to_string(t), int(L / h) / 2);

		timer.Tock(idtpress2);

		timer.Tick(idtuvwmom3);

		fdcolfield.CalcdphidtUmomDiv(fdcolfield.cdUdt2, aa, bb, cc);
		
		fdcolfield.CalcdphidtWmomDiv(fdcolfield.cdWdt2, aa, bb, cc);

		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cUnp(i, j, k) =
						fdcolfield.cUn(i, j, k) +
						(dt)*fdcolfield.cdUdt2(i, j, k);
				}
			}
		}
				
		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cWnp(i, j, k) =
						fdcolfield.cWn(i, j, k) +
						(dt)*fdcolfield.cdWdt2(i, j, k);
				}
			}
		}

		timer.Tock(idtuvwmom3);

		// calculate pressure and enforce continuity

		timer.Tick(idtpress3);

		fdcolfield.EnforceContinuityVelnp(dt, aa, bb, cc, ptr, col, val, rhs, true, lssolver);

		timer.Tock(idtpress3);

		timer.Tick(idtuvwmom4);

		fdcolfield.CalcdphidtUmomDiv(fdcolfield.cdUdt3, aa, bb, cc);

		fdcolfield.CalcdphidtWmomDiv(fdcolfield.cdWdt3, aa, bb, cc);


		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cUnp(i, j, k) =
						fdcolfield.cUn(i, j, k) +
						(dt / 6.0)*(
							fdcolfield.cdUdtn(i, j, k) +
							2.0 * fdcolfield.cdUdt1(i, j, k) +
							2.0 * fdcolfield.cdUdt2(i, j, k) +
							fdcolfield.cdUdt3(i, j, k)
							);
				}
			}
		}
		
		timer.Tock(idtuvwmom4);

		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd)
					{
						continue;
					}
					fdcolfield.cWnp(i, j, k) =
						fdcolfield.cWn(i, j, k) +
						(dt / 6.0)*(
							fdcolfield.cdWdtn(i, j, k) +
							2.0 * fdcolfield.cdWdt1(i, j, k) +
							2.0 * fdcolfield.cdWdt2(i, j, k) +
							fdcolfield.cdWdt3(i, j, k)
							);
				}
			}
		}


		// calculate pressure and enforce continuity

		timer.Tick(idtpress4);

		fdcolfield.EnforceContinuityVelnp(dt, aa, bb, cc, ptr, col, val, rhs, true, lssolver);

		timer.Tock(idtpress4);

		// START SAVE TO FILE

		{
			nn = fdcolfield.CreatedphidxStaggeredSystemTDSimplified(fdcolfield.cUnp, fdcolfield.cNXU, fdcolfield.cNZU,
				aa, bb, cc, fdcolfield.cdphidx);

			lssolver.SolveTridiagonalDestructive(fdcolfield.cdphidx, nn, aa, bb, cc);
			
			nn = fdcolfield.CreatedphidzStaggeredSystemTDSimplified(fdcolfield.cWnp, fdcolfield.cNXW, fdcolfield.cNZW,
				aa, bb, cc, fdcolfield.cdphidz);

			lssolver.SolveTridiagonalDestructive(fdcolfield.cdphidz, nn, aa, bb, cc);
		}

		if ((it % saveEvery) == 0)
		{
			fdcolfield.SaveToCSV("flatplateturb/save" + std::to_string(t));
		}

		double maxDiv = 0.0;

		int maxDivi, maxDivj, maxDivk;
		maxDivi = maxDivj = maxDivk = 0;

		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZP - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXP - 2; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd &&
						k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd)
					{
						continue;
					}

					const double div = fdcolfield.cdphidx(i-1, j, k) +
						 +
						fdcolfield.cdphidz(i, j, k -1);

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
		fdcolfield.cWn = fdcolfield.cWnp;

		t += dt;
	}
}