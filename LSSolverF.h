#pragma once

#define AMGCL_NO_BOOST
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#include <amgcl/preconditioner/dummy.hpp>
#include "preconditionerDummy2.hpp"

typedef amgcl::backend::builtin<float> BackendF;

typedef amgcl::make_solver<
	// Use AMG as preconditioner:
	amgcl::amg<
	BackendF,
	amgcl::coarsening::smoothed_aggregation,
	amgcl::relaxation::spai0
	>,
	// And BiCGStab as iterative solver:
	amgcl::solver::bicgstab<BackendF>
> SolverF;

typedef amgcl::make_solver<
	amgcl::preconditioner::dummy2<BackendF>,
	// And BiCGStab as iterative solver:
	amgcl::solver::bicgstab<BackendF>
> SolverFWithoutAMG;

class LSSolverF
{
public:
	LSSolverF();
	~LSSolverF();

	void SolveSparseCRSWithAMG(std::vector<float>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<float>& val, const std::vector<float>& rhs, float& errorOut, int& nItOut, float tol = 1e-9);

	void PrecondtionCRS(int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<float>& val, float tol = 1e-9);

	void SolvePreconditionedCRS(std::vector<float>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<float>& val, const std::vector<float>& rhs, float& errorOut, int& nItOut);

	void SolveSparseCRS(std::vector<float>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
		const std::vector<float>& val, const std::vector<float>& rhs, float& errorOut, int& nItOut, float tol = 1e-9);

	// solve [aa bb cc] x = b; -> b é a saida e cc será modificado
	static void SolveTridiagonalDestructive(std::vector<float>& b, const size_t N, const std::vector<float>& aa,
		const std::vector<float>& bb, std::vector<float>& cc);

	// solve [aai bbi cci] x = bi; -> b é a saida e cc será modificado
	void SolveTridiagonalDestructiveAVX2(const std::vector<std::vector<float>*>& b,
		const size_t N, 
		const std::vector<const std::vector<float>*>& aa,
		const std::vector<const std::vector<float>*>& bb, 
		const std::vector<std::vector<float>*>& cc);

private:
	SolverF* cPSolverFPreconditioned;
};