#pragma once
#include <iostream>
#include <thread>

#include "SquareCylinder2DStaggeredField.h"
#include "LSSolverF.h"

#include "Timer.h"

void SquareRectAdv2D()
{
	constexpr double h = 0.0002; // passo no espa?o
	constexpr double l = 0.02;   // dimens?o do cilindro
	constexpr double L = 26.0*l; // dimens?o horizontal do campo
	constexpr double H = 15.0*l; // dimens?o vertical do campo
	constexpr int NX = int(L / h + 2); // num. pontos horizontal
	constexpr int NZ = int(H / h + 2); // num. pontos vertical
	constexpr int NN = int(NX * NZ);   // total de pontos
	constexpr double ufree = 0.15; // vel. horizontal farfield
	constexpr double wfree = 0.0;  // vel. vertical farfield
	constexpr double dt = h / 2; // passo no tempo
	constexpr double cfl = 4 * ufree * dt / h; // CFL < 0.3
	constexpr double rho = 1.0; // massa especifica
	constexpr double mu = 3.0e-5;// visc. dinamica
	constexpr double tfinal = 15.0; // tempo final
	constexpr int NT = int(tfinal / dt) + 1; // quantidade de passos no tempo

	constexpr double CellRe = rho * h*ufree / mu; 

	constexpr double Re = rho * l*ufree / mu;

	constexpr int saveEvery = 366 * (0.001/ dt);

	constexpr double dtSave = saveEvery * dt;

	constexpr double conditionForStability = 4 * dt / (Re*h*h);

	std::vector<double> aa, bb, cc, aa1, bb1, cc1, aa2, bb2, cc2;
	Field3DIKJ rhs;
	std::vector<double> val;
	std::vector<int> ptr, col;

	LSSolver lssolver;

	SquareCylinder2DStaggeredField fdcolfield(rho, mu, h, l, L, H, ufree);

	fdcolfield.SetBCFlatPlate();

	std::cout << "Reiniciar ? (y/n)" << std::endl;
	char c = 0;
	std::cin >> c;


	double t = 0.0;

	if (c == 'y')
	{
		std::string basefilename = "";

		std::cout << "Qual o nome do ultimo arquivo? (inclua tudo menos U/V/W/P e .dat)" << std::endl;

		std::cin >> basefilename;

		fdcolfield.ReadField(basefilename);

		std::cout << "Qual o tempo atual da simulacao?" << std::endl;

		std::cin >> t;
	}

	fdcolfield.SaveToCSV("BC");

	std::cout << "cNXCy " << fdcolfield.cNXCy << std::endl;
	std::cout << "cNZCy " << fdcolfield.cNZCy << std::endl;

	std::cout << "lx " << fdcolfield.cNXCy * fdcolfield.ch << std::endl;
	std::cout << "lz " << fdcolfield.cNZCy * fdcolfield.ch << std::endl;

	std::cout << "NNCy " << fdcolfield.cNXCy*fdcolfield.cNZCy << std::endl;

	std::cout << "Re " << rho * ufree * fdcolfield.cNZCy * fdcolfield.ch / mu << std::endl;

	Timer timer;

	const int idttimestep = timer.SetTimer("Complete timestep");

	const int idtuvwmom1 = timer.SetTimer("uvwmom 1 calculation + sum");
	const int idtuvwmom1c = timer.SetTimer("uvwmom 1 sum");
	const int idtpress1 = timer.SetTimer("pressure 1 calculation and continuity");

	const int idtuvwmom2 = timer.SetTimer("uvwmom 2 calculation");
	const int idtuvwmom2c = timer.SetTimer("uvwmom 2 sum");
	const int idtpress2 = timer.SetTimer("pressure 2 calculation and continuity");

	const int idtuvwmom3 = timer.SetTimer("uvwmom 3 calculation");
	const int idtuvwmom3c = timer.SetTimer("uvwmom 3 sum");
	const int idtpress3 = timer.SetTimer("pressure 3 calculation and continuity");

	const int idtuvwmom4 = timer.SetTimer("uvwmom 4 calculation");
	const int idtuvwmom4c = timer.SetTimer("uvwmom 4 sum");
	const int idtpress4 = timer.SetTimer("pressure 4 calculation and continuity");

	const int idtumom = timer.SetTimer("u mom calc");
	const int idtvmom = timer.SetTimer("v mom calc");
	const int idtwmom = timer.SetTimer("w mom calc");

	bool show = true;

	for (int it = 0; it < NT; it++)
	{
		int nn = 0;

		// Set inlet pertubation
		if (t < 10.0)
			for (int j = 0; j < 1; j++)
			{
				for (int k = 1; k < fdcolfield.cNZU - 1; k++)
				{
					fdcolfield.cUnp(0, j, k) =
						fdcolfield.cUfreestream + (sin(t + (k * fdcolfield.ch)) * fdcolfield.cUfreestream * 0.01);
				}
			}
		else
			for (int j = 0; j < 1; j++)
			{
				for (int k = 1; k < fdcolfield.cNZU - 1; k++)
				{
					fdcolfield.cUnp(0, j, k) =
						fdcolfield.cUfreestream + (sin(t + (k * fdcolfield.ch)) * fdcolfield.cUfreestream * 0.00);
				}
			}

		fdcolfield.UpdateBCFlatPlate();

		timer.Tick(idttimestep);

		// momentum eqts.

		timer.Tick(idtuvwmom1);

		timer.Tick(idtumom);

		fdcolfield.InterpolateUWInRespDirections(aa, bb, cc, aa1, bb1, cc1, aa2, bb2, cc2);

		/*fdcolfield.CalcdphidtUmom(fdcolfield.cdUdtn, aa, bb, cc);*/

		std::thread dphidtumomThread(&SquareCylinder2DStaggeredField::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdtn),
			std::ref(aa), std::ref(bb), std::ref(cc));

		timer.Tock(idtumom);

		timer.Tick(idtvmom);

		//fdcolfield.CalcdphidtVmom(fdcolfield.cdVdtn, lssolver, aa1, bb1, cc1);
		
		timer.Tock(idtvmom);

		timer.Tick(idtwmom);

		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdtn, lssolver, aa2, bb2, cc2);

		std::thread dphidtwmomThread(&SquareCylinder2DStaggeredField::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdtn),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));

		timer.Tock(idtwmom);

		dphidtumomThread.join();
		dphidtwmomThread.join();

		timer.Tick(idtuvwmom1c);

#pragma omp parallel for
		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
						)
					{
						continue;
					}

					fdcolfield.cUnp(i, j, k) =
						fdcolfield.cUn(i, j, k) +
						(dt / 2.0)*fdcolfield.cdUdtn(i, j, k);
				}
			}
		}

#pragma omp parallel for
		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
			{
				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
				{
					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd
						)
					{
						continue;
					}
					fdcolfield.cWnp(i, j, k) =
						fdcolfield.cWn(i, j, k) +
						(dt / 2.0)*fdcolfield.cdWdtn(i, j, k);
				}
			}
		}

		timer.Tock(idtuvwmom1c);


		timer.Tock(idtuvwmom1);

		// calculate pressure and enforce continuity

		//fdcolfield.SaveIKCutToCSV("flatplateturb/prepress" + std::to_string(t), int(W / h) / 2);

		timer.Tick(idtpress1);

		//fdcolfield.EnforceContinuityVelnp(dt / 2.0, aa, bb, cc, ptr, col, val, rhs, show, lssolver);

		timer.Tock(idtpress1);

		timer.Tick(idtuvwmom2);
		
		fdcolfield.InterpolateUWInRespDirections(aa, bb, cc, aa1, bb1, cc1, aa2, bb2, cc2);
		
		//fdcolfield.SaveIKCutToCSV("flatplateturb/pospress" + std::to_string(t), int(W / h) / 2);
		
		//fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt1, aa, bb, cc);
		
		std::thread dphidtumomThread1(&SquareCylinder2DStaggeredField::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdt1),
			std::ref(aa), std::ref(bb), std::ref(cc));
		
		
		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt1, aa, bb, cc);
		
		std::thread dphidtwmomThread1(&SquareCylinder2DStaggeredField::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdt1),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));
		
		dphidtumomThread1.join();
		dphidtwmomThread1.join();
		
		timer.Tick(idtuvwmom2c);
		

		for (int j = 0; j < 1; j++)
		{
#pragma omp parallel for
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
#pragma omp parallel for
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
		
		timer.Tock(idtuvwmom2c);
		
		timer.Tock(idtuvwmom2);
		
		//fdcolfield.SaveIKCutToCSV("flatplateturb/prepress" + std::to_string(t), int(W / h) / 2);
		
		// calculate pressure and enforce continuity
		
		timer.Tick(idtpress2);
		
				
		//fdcolfield.EnforceContinuityVelnp(dt / 2.0, aa, bb, cc, ptr, col, val, rhs, show, lssolver);
		
		//fdcolfield.SaveIKCutToCSV("flatplateturb/pospress" + std::to_string(t), int(W / h) / 2);
		
		timer.Tock(idtpress2);
		
		timer.Tick(idtuvwmom3);
		
		fdcolfield.InterpolateUWInRespDirections(aa, bb, cc, aa1, bb1, cc1, aa2, bb2, cc2);
		
		//fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt2, aa, bb, cc);
		
		std::thread dphidtumomThread2(&SquareCylinder2DStaggeredField::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdt2),
			std::ref(aa), std::ref(bb), std::ref(cc));
				
		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt2, aa, bb, cc);
		
		std::thread dphidtwmomThread2(&SquareCylinder2DStaggeredField::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdt2),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));
		
		dphidtumomThread2.join();
		dphidtwmomThread2.join();
		
		timer.Tick(idtuvwmom3c);
		
		for (int j = 0; j < 1; j++)
		{

#pragma omp parallel for
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
#pragma omp parallel for
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
		
		timer.Tock(idtuvwmom3c);
		
		timer.Tock(idtuvwmom3);
		
		// calculate pressure and enforce continuity
		
		timer.Tick(idtpress3);
		
		fdcolfield.InterpolateUWInRespDirections(aa, bb, cc, aa1, bb1, cc1, aa2, bb2, cc2);
		
		//fdcolfield.SaveIKCutToCSV("flatplateturb/prepress2" + std::to_string(t), int(W / h) / 2);
		
		//fdcolfield.EnforceContinuityVelnp(dt, aa, bb, cc, ptr, col, val, rhs, show, lssolver);
		
		//fdcolfield.SaveIKCutToCSV("flatplateturb/pospress2" + std::to_string(t), int(W / h) / 2);
		
		timer.Tock(idtpress3);
		
		timer.Tick(idtuvwmom4);
		
		//fdcolfield.CalcdphidtUmom(fdcolfield.cdUdt3, aa, bb, cc);
		
		std::thread dphidtumomThread3(&SquareCylinder2DStaggeredField::CalcdphidtUmomAdv, &fdcolfield, std::ref(fdcolfield.cdUdt3),
			std::ref(aa), std::ref(bb), std::ref(cc));
				
		//fdcolfield.CalcdphidtWmom(fdcolfield.cdWdt3, aa, bb, cc);
		
		std::thread dphidtwmomThread3(&SquareCylinder2DStaggeredField::CalcdphidtWmomAdv, &fdcolfield, std::ref(fdcolfield.cdWdt3),
			std::ref(aa2), std::ref(bb2), std::ref(cc2));
		
		dphidtumomThread3.join();
		dphidtwmomThread3.join();
		
		timer.Tick(idtuvwmom4c);
		
//#pragma omp parallel for
//		for (int j = 1; j < fdcolfield.cNYU - 1; j++)
//		{
//			for (int k = 1; k < fdcolfield.cNZU - 1; k++)
//			{
//				for (int i = 1; i < fdcolfield.cNXU - 1; i++)
//				{
//					if (i >= fdcolfield.ciCyInit - 1 && i < fdcolfield.ciCyEnd
//						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
//						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
//					{
//						continue;
//					}
//					fdcolfield.cUnp(i, j, k) =
//						fdcolfield.cUn(i, j, k) +
//						(dt)*fdcolfield.cdUdt3(i, j, k);
//				}
//			}
//		}
//
//#pragma omp parallel for
//		for (int j = 1; j < fdcolfield.cNYV - 1; j++)
//		{
//			for (int k = 1; k < fdcolfield.cNZV - 1; k++)
//			{
//				for (int i = 1; i < fdcolfield.cNXV - 1; i++)
//				{
//					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
//						&& k >= fdcolfield.ckCyInit && k < fdcolfield.ckCyEnd
//						&& j >= fdcolfield.cjCyInit - 1 && j < fdcolfield.cjCyEnd)
//					{
//						continue;
//					}
//					fdcolfield.cVnp(i, j, k) =
//						fdcolfield.cVn(i, j, k) +
//						(dt)*fdcolfield.cdVdt3(i, j, k);
//				}
//			}
//		}
//
//#pragma omp parallel for
//		for (int j = 1; j < fdcolfield.cNYW - 1; j++)
//		{
//			for (int k = 1; k < fdcolfield.cNZW - 1; k++)
//			{
//				for (int i = 1; i < fdcolfield.cNXW - 1; i++)
//				{
//					if (i >= fdcolfield.ciCyInit && i < fdcolfield.ciCyEnd
//						&& k >= fdcolfield.ckCyInit - 1 && k < fdcolfield.ckCyEnd
//						&& j >= fdcolfield.cjCyInit && j < fdcolfield.cjCyEnd)
//					{
//						continue;
//					}
//					fdcolfield.cWnp(i, j, k) =
//						fdcolfield.cWn(i, j, k) +
//						(dt)*fdcolfield.cdWdt3(i, j, k);
//				}
//			}
//		}
		
		
		/*fdcolfield.SaveIKCutToCSV("dudtn", fdcolfield.cdUdtn, fdcolfield.cdUdtn.cNX, fdcolfield.cdUdtn.cNY, fdcolfield.cdUdtn.cNZ,
			fdcolfield.cdUdtn.cNY / 2);
		
		fdcolfield.SaveIKCutToCSV("dudt1", fdcolfield.cdUdt1, fdcolfield.cdUdtn.cNX, fdcolfield.cdUdtn.cNY, fdcolfield.cdUdtn.cNZ,
			fdcolfield.cdUdtn.cNY / 2);
		
		fdcolfield.SaveIKCutToCSV("dudt2", fdcolfield.cdUdt2, fdcolfield.cdUdtn.cNX, fdcolfield.cdUdtn.cNY, fdcolfield.cdUdtn.cNZ,
			fdcolfield.cdUdtn.cNY / 2);
		
		fdcolfield.SaveIKCutToCSV("dudt3", fdcolfield.cdUdt3, fdcolfield.cdUdtn.cNX, fdcolfield.cdUdtn.cNY, fdcolfield.cdUdtn.cNZ,
			fdcolfield.cdUdtn.cNY / 2);*/
		

		for (int j = 0; j < 1; j++)
		{
#pragma omp parallel for
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
#pragma omp parallel for
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
		
		timer.Tock(idtuvwmom4c);

				// calculate pressure and enforce continuity

		timer.Tick(idtpress4);

		//fdcolfield.SaveIKCutToCSV("flatplateturb/prepress3" + std::to_string(t), int(W / h) / 2);

		fdcolfield.EnforceContinuityVelnp(dt, aa, bb, cc, ptr, col, val, rhs, show, lssolver);

		//fdcolfield.SaveIKCutToCSV("flatplateturb/pospress3" + std::to_string(t), int(W / h) / 2);

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

		double maxCFL = 0.0;
		int maxCFLi, maxCFLj, maxCFLk;
		maxCFLi = maxCFLj = maxCFLk = 0;

		for (int j = 0; j < 1; j++)
		{
			for (int k = 1; k < fdcolfield.cWnpc.cNZ; k++)
			{
				for (int i = 1; i < fdcolfield.cUnpc.cNX; i++)
				{
					const double vels = std::fabs(fdcolfield.cUnpc(i - 1, j, k)) +
						std::fabs(fdcolfield.cWnpc(i, j, k - 1));

					const double cfll = vels * dt / h;

					if (maxCFL < cfll)
					{
						maxCFL = cfll;

						maxCFLi = i;
						maxCFLj = j;
						maxCFLk = k;
					}
				}
			}
		}

		if ((it % saveEvery) == 0 || (t > 12.0 && (it % (saveEvery / 6)) == 0))
		{
			fdcolfield.SaveToCSV("flatplateturb/save" + std::to_string(t));

			//fdcolfield.SaveIKCutToCSV("flatplateturb/save" + std::to_string(t), 10);
			//
			//fdcolfield.SaveIKCutToCSV("flatplateturb/save" + std::to_string(t), 30);

			fdcolfield.SaveField("output/field" + std::to_string(t));
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

					const double div = fdcolfield.cdphidx(i - 1, j, k) +
						fdcolfield.cdphidz(i, j, k - 1);

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

		timer.Tock(idttimestep);

		std::cout << "#######################" << std::endl;

		std::cout << "max div" << maxDiv << " i " << maxDivi << " j " << maxDivj << " k " << maxDivk << std::endl;
		std::cout << "max CFL " << maxCFL << " i " << maxCFLi << " j " << maxCFLj << " k " << maxCFLk << std::endl;
		std::cout << "tempo " << t << std::endl;
		std::cout << "tempo ultimo timestep " << timer.GetLastTickTock(idttimestep) << std::endl;
		std::cout << "tempo restante " << (NT - it) * timer.GetLastTickTock(idttimestep) / 60.0 / 60.0 / 24.0 << " dias" << std::endl;

		if (show)
			timer.WriteToCoutAllTickTock();

		std::cout << "#######################" << std::endl;

		if (GetAsyncKeyState('h') || GetAsyncKeyState('H'))
		{
			show = false;

			std::cout << "Hiding timer couts. (H to hide, S to show)" << std::endl;
		}

		if (GetAsyncKeyState('s') || GetAsyncKeyState('S'))
		{
			show = true;

			std::cout << "Showing timer couts.(H to hide, S to show)" << std::endl;
		}

		// END SAVE TO FILE

		// set old to new
		fdcolfield.cUn = fdcolfield.cUnp;
		fdcolfield.cWn = fdcolfield.cWnp;

		t += double(dt);
	}
}