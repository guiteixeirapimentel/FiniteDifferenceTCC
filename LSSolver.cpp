#include "LSSolver.h"
#include <immintrin.h>

void LSSolver::SolveSparseCRS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, double tol)
{
	SolverWithoutAMG::params prm;
	prm.solver.tol = tol;

	SolverWithoutAMG solve(std::tie(nCol, ptr, col, val), prm);
	std::tie(nItOut, errorOut) = solve(std::tie(nCol, ptr, col, val), rhs, std::move(x));
}

LSSolver::LSSolver()
	:
	cPSolverFPreconditioned(nullptr)
{}
LSSolver::~LSSolver()
{
	if (cPSolverFPreconditioned)
	{
		delete cPSolverFPreconditioned;
		cPSolverFPreconditioned = nullptr;
	}
}

void LSSolver::PrecondtionCRS(int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, double tol)
{
	if (cPSolverFPreconditioned)
	{
		delete cPSolverFPreconditioned;
		cPSolverFPreconditioned = nullptr;
	}

	Solver::params prm;
	prm.solver.tol = tol;
	prm.solver.maxiter = 25;
	cPSolverFPreconditioned = new Solver(std::tie(nCol, ptr, col, val), prm);
}

void LSSolver::SolvePreconditionedCRS(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut)
{
	std::tie(nItOut, errorOut) = (*cPSolverFPreconditioned)(std::tie(nCol, ptr, col, val), rhs, std::move(x));
}

void LSSolver::SolveSparseCRSWithAMG(std::vector<double>&& x, int nCol, const std::vector<int>& ptr, std::vector<int>& col,
	const std::vector<double>& val, const std::vector<double>& rhs, double& errorOut, int& nItOut, double tol)
{
	Solver::params prm;
	prm.solver.tol = tol;

	Solver solve(std::tie(nCol, ptr, col, val), prm);

	std::tie(nItOut, errorOut) = solve(rhs, std::move(x));
}

void LSSolver::SolveTridiagonalDestructiveAVX2(const std::vector<std::vector<double>*>& x,
	const size_t N,
	const std::vector <const std::vector<double>*>& aa,
	const std::vector<const std::vector<double>*>& bb,
	const std::vector<std::vector<double>*>& cc)
{
	size_t n;

	/*
	resolve Ax = v onde A é uma matriz tridiagonal composta pelos veores a, b, c
	note que o conteúdo do vetor de entrada c será modificado, tornando esta uma função de um único tempo de uso
	x[] - inicialmente contém o vector de entrada v e retorna a solução x. indexados por[0, ..., N - 1]
	N - número de equações
	a[] - subdiagonal (diagonal abaixo da diagonal principal) -- indexados por [1, ..., N - 1]
	b[] - matriz principal, indexados por [0, ..., N - 1]
	c[] - superdiagonal (diagonal acima da diagonal principal) -- indexedos por [0, ..., N - 2]
	*/

	// 64 bit "double" registers
	__m256 _c, _x, _b, _a, _m, _tmp1, _one;

	_one = _mm256_set1_ps(1.0f);

	// _mm256_div_ps() // divisão de singles 
	// _mm256_mul_ps() // multiplicação de singles
	// _mm256_add_ps() // soma de singles
	// _mm256_sub_ps() // subtração de singles


	//c[0] = c[0] / b[0];

	_c.m256_f32[0] = (*cc[0])[0];
	_c.m256_f32[1] = (*cc[1])[0];
	_c.m256_f32[2] = (*cc[2])[0];
	_c.m256_f32[3] = (*cc[3])[0];
	_c.m256_f32[4] = (*cc[4])[0];
	_c.m256_f32[5] = (*cc[5])[0];
	_c.m256_f32[6] = (*cc[6])[0];
	_c.m256_f32[7] = (*cc[7])[0];

	_b.m256_f32[0] = (*bb[0])[0];
	_b.m256_f32[1] = (*bb[1])[0];
	_b.m256_f32[2] = (*bb[2])[0];
	_b.m256_f32[3] = (*bb[3])[0];
	_b.m256_f32[4] = (*bb[4])[0];
	_b.m256_f32[5] = (*bb[5])[0];
	_b.m256_f32[6] = (*bb[6])[0];
	_b.m256_f32[7] = (*bb[7])[0];

	_c = _mm256_div_ps(_c, _b);

	//x[0] = x[0] / b[0];

	_x.m256_f32[0] = (*cc[0])[0];
	_x.m256_f32[1] = (*cc[1])[0];
	_x.m256_f32[2] = (*cc[2])[0];
	_x.m256_f32[3] = (*cc[3])[0];
	_x.m256_f32[4] = (*cc[4])[0];
	_x.m256_f32[5] = (*cc[5])[0];
	_x.m256_f32[6] = (*cc[6])[0];
	_x.m256_f32[7] = (*cc[7])[0];

	//_b.m256_f32[0] = bb[0][0];
	//_b.m256_f32[1] = bb[1][0];
	//_b.m256_f32[2] = bb[2][0];
	//_b.m256_f32[3] = bb[3][0];
	//_b.m256_f32[4] = bb[4][0];
	//_b.m256_f32[5] = bb[5][0];
	//_b.m256_f32[6] = bb[6][0];
	//_b.m256_f32[7] = bb[7][0];

	_x = _mm256_div_ps(_x, _b);

	(*x[0])[0] = _x.m256_f32[0];
	(*x[1])[0] = _x.m256_f32[1];
	(*x[2])[0] = _x.m256_f32[2];
	(*x[3])[0] = _x.m256_f32[3];
	(*x[4])[0] = _x.m256_f32[4];
	(*x[5])[0] = _x.m256_f32[5];
	(*x[6])[0] = _x.m256_f32[6];
	(*x[7])[0] = _x.m256_f32[7];

	/* loop de 1 a N - 1 inclusive */
	for (n = 1; n < N; n++)
	{
		//double m = 1.0 / (b[n] - a[n] * c[n - 1]);

		// m = a[n] * c[n - 1];
		_a.m256_f32[0] = (*aa[0])[n];
		_a.m256_f32[1] = (*aa[1])[n];
		_a.m256_f32[2] = (*aa[2])[n];
		_a.m256_f32[3] = (*aa[3])[n];
		_a.m256_f32[4] = (*aa[4])[n];
		_a.m256_f32[5] = (*aa[5])[n];
		_a.m256_f32[6] = (*aa[6])[n];
		_a.m256_f32[7] = (*aa[7])[n];

		_c.m256_f32[0] = (*cc[0])[n - 1];
		_c.m256_f32[1] = (*cc[1])[n - 1];
		_c.m256_f32[2] = (*cc[2])[n - 1];
		_c.m256_f32[3] = (*cc[3])[n - 1];
		_c.m256_f32[4] = (*cc[4])[n - 1];
		_c.m256_f32[5] = (*cc[5])[n - 1];
		_c.m256_f32[6] = (*cc[6])[n - 1];
		_c.m256_f32[7] = (*cc[7])[n - 1];

		_m = _mm256_mul_ps(_a, _c);

		// m = b[n] - (a[n] * c[n - 1])
		_b.m256_f32[0] = (*bb[0])[n];
		_b.m256_f32[1] = (*bb[1])[n];
		_b.m256_f32[2] = (*bb[2])[n];
		_b.m256_f32[3] = (*bb[3])[n];
		_b.m256_f32[4] = (*bb[4])[n];
		_b.m256_f32[5] = (*bb[5])[n];
		_b.m256_f32[6] = (*bb[6])[n];
		_b.m256_f32[7] = (*bb[7])[n];
		_m = _mm256_sub_ps(_b, _m);

		// m = 1 / (b[n] - a[n] * c[n - 1])
		_m = _mm256_div_ps(_one, _m);


		//c[n] = c[n] * m;

		_c.m256_f32[0] = (*cc[0])[n];
		_c.m256_f32[1] = (*cc[1])[n];
		_c.m256_f32[2] = (*cc[2])[n];
		_c.m256_f32[3] = (*cc[3])[n];
		_c.m256_f32[4] = (*cc[4])[n];
		_c.m256_f32[5] = (*cc[5])[n];
		_c.m256_f32[6] = (*cc[6])[n];
		_c.m256_f32[7] = (*cc[7])[n];

		_c = _mm256_mul_ps(_c, _m);

		//x[n] = (x[n] - a[n] * x[n - 1]) * m;

		_x.m256_f32[0] = (*x[0])[n - 1];
		_x.m256_f32[1] = (*x[1])[n - 1];
		_x.m256_f32[2] = (*x[2])[n - 1];
		_x.m256_f32[3] = (*x[3])[n - 1];
		_x.m256_f32[4] = (*x[4])[n - 1];
		_x.m256_f32[5] = (*x[5])[n - 1];
		_x.m256_f32[6] = (*x[6])[n - 1];
		_x.m256_f32[7] = (*x[7])[n - 1];

		_tmp1 = _mm256_mul_ps(_a, _x);

		_x.m256_f32[0] = (*x[0])[n];
		_x.m256_f32[1] = (*x[1])[n];
		_x.m256_f32[2] = (*x[2])[n];
		_x.m256_f32[3] = (*x[3])[n];
		_x.m256_f32[4] = (*x[4])[n];
		_x.m256_f32[5] = (*x[5])[n];
		_x.m256_f32[6] = (*x[6])[n];
		_x.m256_f32[7] = (*x[7])[n];

		_x = _mm256_sub_ps(_x, _tmp1);

		_x = _mm256_mul_ps(_x, _m);

		(*x[0])[n] = _x.m256_f32[0];
		(*x[1])[n] = _x.m256_f32[1];
		(*x[2])[n] = _x.m256_f32[2];
		(*x[3])[n] = _x.m256_f32[3];
		(*x[4])[n] = _x.m256_f32[4];
		(*x[5])[n] = _x.m256_f32[5];
		(*x[6])[n] = _x.m256_f32[6];
		(*x[7])[n] = _x.m256_f32[7];
	}

	/* loop de N - 2 a 0 inclusive */
	for (n = N - 1; n-- > 0; )
	{
		//x[n] = x[n] - c[n] * x[n + 1];
		_x.m256_f32[0] = (*x[0])[n + 1];
		_x.m256_f32[1] = (*x[1])[n + 1];
		_x.m256_f32[2] = (*x[2])[n + 1];
		_x.m256_f32[3] = (*x[3])[n + 1];
		_x.m256_f32[4] = (*x[4])[n + 1];
		_x.m256_f32[5] = (*x[5])[n + 1];
		_x.m256_f32[6] = (*x[6])[n + 1];
		_x.m256_f32[7] = (*x[7])[n + 1];
		_tmp1 = _mm256_mul_ps(_c, _x);

		_x.m256_f32[0] = (*x[0])[n];
		_x.m256_f32[1] = (*x[1])[n];
		_x.m256_f32[2] = (*x[2])[n];
		_x.m256_f32[3] = (*x[3])[n];
		_x.m256_f32[4] = (*x[4])[n];
		_x.m256_f32[5] = (*x[5])[n];
		_x.m256_f32[6] = (*x[6])[n];
		_x.m256_f32[7] = (*x[7])[n];

		_x = _mm256_sub_ps(_x, _tmp1);

		(*x[0])[n] = _x.m256_f32[0];
		(*x[1])[n] = _x.m256_f32[1];
		(*x[2])[n] = _x.m256_f32[2];
		(*x[3])[n] = _x.m256_f32[3];
		(*x[4])[n] = _x.m256_f32[4];
		(*x[5])[n] = _x.m256_f32[5];
		(*x[6])[n] = _x.m256_f32[6];
		(*x[7])[n] = _x.m256_f32[7];
	}
}

void LSSolver::SolveTridiagonalDestructive(std::vector<double>& x, const size_t N, const std::vector<double>& a,
	const std::vector<double>& b, std::vector<double>& c)
{
	size_t n;

	/*
	resolve Ax = v onde A é uma matriz tridiagonal composta pelos veores a, b, c
	note que o conteúdo do vetor de entrada c será modificado, tornando esta uma função de um único tempo de uso
	x[] - inicialmente contém o vector de entrada v e retorna a solução x. indexados por[0, ..., N - 1]
	N - número de equações
	a[] - subdiagonal (diagonal abaixo da diagonal principal) -- indexados por [1, ..., N - 1]
	b[] - matriz principal, indexados por [0, ..., N - 1]
	c[] - superdiagonal (diagonal acima da diagonal principal) -- indexedos por [0, ..., N - 2]
	*/

	c[0] = c[0] / b[0];
	x[0] = x[0] / b[0];

	/* loop de 1 a N - 1 inclusive */
	for (n = 1; n < N; n++) {
		double m = 1.0 / (b[n] - a[n] * c[n - 1]);
		c[n] = c[n] * m;
		x[n] = (x[n] - a[n] * x[n - 1]) * m;
	}

	/* loop de N - 2 a 0 inclusive */
	for (n = N - 1; n-- > 0; )
		x[n] = x[n] - c[n] * x[n + 1];
}

void LSSolver::SolveTridiagonalDestructive(Field3D& x, const size_t N, const std::vector<double>& a,
	const std::vector<double>& b, std::vector<double>& c)
{
	size_t n;

	/*
	resolve Ax = v onde A é uma matriz tridiagonal composta pelos veores a, b, c
	note que o conteúdo do vetor de entrada c será modificado, tornando esta uma função de um único tempo de uso
	x[] - inicialmente contém o vector de entrada v e retorna a solução x. indexados por[0, ..., N - 1]
	N - número de equações
	a[] - subdiagonal (diagonal abaixo da diagonal principal) -- indexados por [1, ..., N - 1]
	b[] - matriz principal, indexados por [0, ..., N - 1]
	c[] - superdiagonal (diagonal acima da diagonal principal) -- indexedos por [0, ..., N - 2]
	*/

	c[0] = c[0] / b[0];
	x[0] = x[0] / b[0];

	/* loop de 1 a N - 1 inclusive */
	for (n = 1; n < N; n++) {
		double m = 1.0 / (b[n] - a[n] * c[n - 1]);
		c[n] = c[n] * m;
		x[n] = (x[n] - a[n] * x[n - 1]) * m;
	}

	/* loop de N - 2 a 0 inclusive */
	for (n = N - 1; n-- > 0; )
		x[n] = x[n] - c[n] * x[n + 1];
}