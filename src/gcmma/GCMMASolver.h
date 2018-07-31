#pragma once

#include <vector>

/*

BETA VERSION  0.99

Sequential MMA solver using a dual interior point method

Code by Niels Aage, February 2013

The class solves a general non-linear programming problem
on standard from, i.e. non-linear objective f0, m non-linear
inequality constraints f and box constraints on the n
design variables xmin, xmax.

       min_x^n f(x)
       s.t. g_j(x) < 0,   j = 1,m
       xmin < x_i < xmax, i = 1,n

Each call to OuterUpdate() sets up and solve the following
convex subproblem:

  min_x     sum(p0j./(U-x)+q0j./(x-L)) + a0*z + sum(c.*y + 0.5*d.*y.^2)

  s.t.      sum(pij./(U-x)+qij./(x-L)) - ai*z - yi <= bi, i = 1,m
            Lj < alphaj <=  xj <= betaj < Uj,  j = 1,n
            yi >= 0, i = 1,m
            z >= 0.

NOTE: a0 == 1 in this implementation !!!!

*/

class GCMMASolver {

public:
	GCMMASolver(int n,int m, double a = 0.0, double c = 1000.0, double d = 0.0);

	void SetAsymptotes(double init, double decrease, double increase);

	// Compute [L, U, raa0, raa], build and solve GCMMA subproblem, compute [f0app, fapp]
	void OuterUpdate(double *xmma, const double *xval, double f0x, const double *df0dx,
		const double *fx, const double *dfdx, const double *xmin, const double *xmax);

	// Update [raa0, raa], build and solve GCMMA subproblem, compute [f0app, fapp]
	void InnerUpdate(double *xmma, double f0xnew, const double *fxnew,
		const double *xval, double f0x, const double *df0dx, const double *fx,
		const double *dfdx, const double *xmin, const double *xmax);

	// Check whether the new solution is conservative
	bool ConCheck(double f0xnew, const double *fxnew) const;

	void Reset() { outeriter = 0; };

private:
	int n, m, outeriter;

	const double raa0eps;
	const double raaeps;
	const double xmamieps;
	const double epsimin;

	const double move, albefa;
	double asyminit, asymdec, asyminc;

	double raa0;
	std::vector<double> raa;

	std::vector<double> a, c, d;
	std::vector<double> y;
	double z;

	std::vector<double> lam, mu, s;
	std::vector<double> low, upp, alpha, beta, p0, q0, pij, qij, b, grad, hess;

	double r0, f0app;
	std::vector<double> r, fapp;

	std::vector<double> xold1, xold2;

private:
	// Compute [low, upp, raa0, raa]
	void Asymp(const double *xval, const double *df0dx,
		const double *dfdx, const double *xmin, const double *xmax);

	// Update [raa0, raa]
	void RaaUpdate(const double *xmma, const double *xval, double f0xnew,
		const double *fxnew, const double *xmin, const double *xmax);

	// Build CGMMA subproblem
	void GenSub(const double *xval, double f0x, const double *df0dx, const double *fx,
		const double *dfdx, const double *xmin, const double *xmax);

	// Compute [f0app, fapp]
	void ComputeApprox(const double *xmma);

	void SolveDSA(double *x);
	void SolveDIP(double *x);

	void XYZofLAMBDA(double *x);

	void DualGrad(double *x);
	void DualHess(double *x);
	void DualLineSearch();
	double DualResidual(double *x, double epsi);

	static void Factorize(double *K, int n);
	static void Solve(double *K, double *x, int n);
};
