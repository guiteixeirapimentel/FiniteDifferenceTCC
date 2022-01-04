#pragma once
#include <iostream>
#include <thread>

#include "SquareCylinder3DStaggeredFieldF.h"
#include "LSSolverF.h"

#include "Timer.h"

void SquareRectFAdv()
{
	constexpr float h = 0.002;

	constexpr float l = 10.0*h;

	constexpr float L = 26.0*l;
	constexpr float H = 12.0*l;
	constexpr float W = 13.0*l;

	constexpr int NX = int(L / h + 2);
	constexpr int NY = int(W / h + 2);
	constexpr int NZ = int(H / h + 2);

	constexpr int NN = int(NX * NY * NZ);

	constexpr float ufree = 0.15;
	constexpr float vfree = 0.0;
	constexpr float wfree = 0.0;

	constexpr float dt = h / 2 * 0.8;

	constexpr float cfl = 4*ufree * dt / h;

	constexpr float rho = 1.0;
	constexpr float mu = 3e-5;

	constexpr float tfinal = 300.0;

	constexpr int NT = int(tfinal / dt) + 1;

	constexpr float CellRe = rho * h*ufree / mu;

	constexpr float Re = rho * l*ufree / mu;

	constexpr int saveEvery = 10;

	constexpr float dtSave = saveEvery * dt;

	constexpr float conditionForStability = 4 * dt / (Re*h*h);

	std::vector<float> aa, bb, cc, aa1, bb1, cc1, aa2, bb2, cc2;
	std::vector<float> rhs;
	std::vector<float> val;
	std::vector<int> ptr, col;

	LSSolverF lssolver;

	SquareCylinder3DStaggeredFieldF fdcolfield(rho, mu, h, l, L, H, W, ufree);

	fdcolfield.SetBCFlatPlate();

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

	const int idtumom = timer.SetTimer("u mom calc");
	const int idtvmom = timer.SetTimer("v mom calc");
	const int idtwmom = timer.SetTimer("w mom calc");

	for (int it = 0; it < NT; it++)
	{
		int nn = 0;

		// Set inlet pertubation
		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				fdcolfield.cUnp[(j * fdcolfield.cNZU * fdcolfield.cNXU) + (k * fdcolfield.cNXU) + 0] =
					fdcolfield.cUfreestream + (sin(t + (k * fdcolfield.ch)) * fdcolfield.cUfreestream * 0.0);
			}
		}

		fdcolfield.UpdateBCFlatPlate();

		timer.Tick(idttimestep);

		// momentum eqts.

		timer.Tick(idtuvwmom1);

		timer.Tick(idtumom);

		/*fdcolfield.CalcdphidtUmom(fdcolfield.cdUdtn, aa, bb, cc);*/

		std::thread dphidtumomThread(&SquareCylinder3DStaggeredFieldF::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdtn),
			std::ref(aa), std::ref(bb), std::ref(cc));

		timer.Tock(idtumom);

		timer.Tick(idtvmom);

		//fdcolfield.CalcdphidtVmom(fdcolfield.cdVdtn, lssolver, aa1, bb1, cc1);

		std::thread dphidtvmomThread(&SquareCylinder3DStaggeredFieldF::CalcdphidtVmomAdv, &fdcolfield, std::ref(fdcolfield.cdVdtn),
			std::ref(aa1), std::ref(bb1), std::ref(cc1));

		timer.Tock(idtvmom);

		timer.Tick(idtwmom);

		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdtn, lssolver, aa2, bb2, cc2);

		std::thread dphidtwmomThread(&SquareCylinder3DStaggeredFieldF::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdtn),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));

		timer.Tock(idtwmom);

		dphidtumomThread.join();
		dphidtvmomThread.join();
		dphidtwmomThread.join();

		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}

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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit - 1 && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
					fdcolfield.cWnp[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] =
						fdcolfield.cWn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
						(dt / 2.0)*fdcolfield.cdWdtn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)];
				}
			}
		}


		timer.Tock(idtuvwmom1);

		// calculate pressure and enforce continuity

		fdcolfield.SaveIKCutToCSV("flatplateturb/prepress" + std::to_string(t), int(W / h) / 2);

		timer.Tick(idtpress1);

		fdcolfield.EnforceContinuityVelnp(dt / 2.0, aa, bb, cc, ptr, col, val, rhs, lssolver);

		timer.Tock(idtpress1);

		timer.Tick(idtuvwmom2);

		fdcolfield.SaveIKCutToCSV("flatplateturb/pospress" + std::to_string(t), int(W / h) / 2);

		//fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt1, aa, bb, cc);

		std::thread dphidtumomThread1(&SquareCylinder3DStaggeredFieldF::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdt1),
			std::ref(aa), std::ref(bb), std::ref(cc));

		//fdcolfield.CalcdphidtVmom(fdcolfield.cdVdt1, aa, bb, cc);

		std::thread dphidtvmomThread1(&SquareCylinder3DStaggeredFieldF::CalcdphidtVmomAdv, &fdcolfield, std::ref(fdcolfield.cdVdt1),
			std::ref(aa1), std::ref(bb1), std::ref(cc1));

		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt1, aa, bb, cc);

		std::thread dphidtwmomThread1(&SquareCylinder3DStaggeredFieldF::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdt1),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));

		dphidtumomThread1.join();
		dphidtvmomThread1.join();
		dphidtwmomThread1.join();

		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit - 1 && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
					fdcolfield.cWnp[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] =
						fdcolfield.cWn[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)] +
						(dt / 2.0)*fdcolfield.cdWdt1[i + (k*fdcolfield.cNXW) + (j*fdcolfield.cNXW*fdcolfield.cNZW)];
				}
			}
		}

		timer.Tock(idtuvwmom2);


		//fdcolfield.SaveIKCutToCSV("flatplateturb/prepress" + std::to_string(t), int(W / h) / 2);

		// calculate pressure and enforce continuity

		timer.Tick(idtpress2);


		fdcolfield.EnforceContinuityVelnp(dt / 2.0, aa, bb, cc, ptr, col, val, rhs, lssolver);

		//fdcolfield.SaveJKCutToCSV("flatplateturb/premom" + std::to_string(t), int(L / h) / 2);

		timer.Tock(idtpress2);

		timer.Tick(idtuvwmom3);

		//fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt2, aa, bb, cc);

		std::thread dphidtumomThread2(&SquareCylinder3DStaggeredFieldF::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdt2),
			std::ref(aa), std::ref(bb), std::ref(cc));

		//fdcolfield.CalcdphidtVmom(fdcolfield.cdVdt2, aa, bb, cc);

		std::thread dphidtvmomThread2(&SquareCylinder3DStaggeredFieldF::CalcdphidtVmomAdv, &fdcolfield, std::ref(fdcolfield.cdVdt2),
			std::ref(aa1), std::ref(bb1), std::ref(cc1));

		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt2, aa, bb, cc);

		std::thread dphidtwmomThread2(&SquareCylinder3DStaggeredFieldF::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdt2),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));

		dphidtumomThread2.join();
		dphidtvmomThread2.join();
		dphidtwmomThread2.join();

		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit - 1 && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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

		//fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt3, aa, bb, cc);

		std::thread dphidtumomThread3(&SquareCylinder3DStaggeredFieldF::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdt3),
			std::ref(aa), std::ref(bb), std::ref(cc));

		//fdcolfield.CalcdphidtVmom(fdcolfield.cdVdt3, aa, bb, cc);

		std::thread dphidtvmomThread3(&SquareCylinder3DStaggeredFieldF::CalcdphidtVmomAdv, &fdcolfield, std::ref(fdcolfield.cdVdt3),
			std::ref(aa1), std::ref(bb1), std::ref(cc1));

		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt3, aa, bb, cc);

		std::thread dphidtwmomThread3(&SquareCylinder3DStaggeredFieldF::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdt3),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));

		dphidtumomThread3.join();
		dphidtvmomThread3.join();
		dphidtwmomThread3.join();

		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit - 1 && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd
						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
					{
						continue;
					}
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
			fdcolfield.SaveIKCutToCSV("flatplateturb/save" + std::to_string(t), fdcolfield.cNYP / 2);

			//fdcolfield.SaveIKCutToCSV("flatplateturb/save" + std::to_string(t), 10);
			//
			//fdcolfield.SaveIKCutToCSV("flatplateturb/save" + std::to_string(t), 30);

			fdcolfield.SaveJKCutToCSV("flatplateturb/save" + std::to_string(t), fdcolfield.cNXCy * 3 + 2);//fdcolfield.cNXP/2);

			fdcolfield.SaveIJCutToCSV("flatplateturb/save" + std::to_string(t), fdcolfield.cNZP / 2);

			fdcolfield.SaveField("data/field" + std::to_string(t));
		}

		float maxDiv = 0.0;

		int maxDivi, maxDivj, maxDivk;
		maxDivi = maxDivj = maxDivk = 0;

		for (int j = 1; j < fdcolfield.cNYP - 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZP - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXP - 2; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd && j >= fdcolfield.cjCyInit &&
						j < fdcolfield.cjCyEnd && k >= fdcolfield.ckCyInit && k < fdcolfield.cjCyEnd)
					{
						continue;
					}

					const float div = fdcolfield.cdphidx[i - 1 + ((k)*(fdcolfield.cNXU - 1)) + ((j) * (fdcolfield.cNXU - 1) * fdcolfield.cNZU)] +
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

		t += double(dt);
	}
}