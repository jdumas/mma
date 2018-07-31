#pragma once

#include <vector>

/*

BETA VERSION  0.99

Sequential MMA solver using a dual interior point method

Code by Niels Aage, February 2013

The class solves a general non-linear programming problem
on standard from, i.e. non-linear objective f, m non-linear
inequality constraints g and box constraints on the n
design variables xmin, xmax.

       min_x^n f(x)
       s.t. g_j(x) < 0,   j = 1,m
       xmin < x_i < xmax, i = 1,n

Each call to Update() sets up and solve the following
convex subproblem:

  min_x     sum(p0j./(U-x)+q0j./(x-L)) + a0*z + sum(c.*y + 0.5*d.*y.^2)

  s.t.      sum(pij./(U-x)+qij./(x-L)) - ai*z - yi <= bi, i = 1,m
            Lj < alphaj <=  xj <= betaj < Uj,  j = 1,n
            yi >= 0, i = 1,m
            z >= 0.

NOTE: a0 == 1 in this implementation !!!!

*/

class MMASolver {

  public:
	MMASolver(int n, int m, double a = 0.0, double c = 1000.0, double d = 0.0);

	void SetAsymptotes(double init, double decrease, double increase);

	void ConstraintModification(bool conMod) {}

	void Update(double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
	            const double *xmax);

	void Reset() { iter = 0; };

  private:
	int n, m, iter;

	const double xmamieps;
	const double epsimin;

	const double raa0;
	const double move, albefa;
	double asyminit, asymdec, asyminc;

	std::vector<double> a, c, d;
	std::vector<double> y;
	double z;

	std::vector<double> lam, mu, s;
	std::vector<double> low, upp, alpha, beta, p0, q0, pij, qij, b, grad, hess;

	std::vector<double> xold1, xold2;

	void GenSub(const double *xval, const double *dfdx, const double *gx, const double *dgdx, const double *xmin,
	            const double *xmax);

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
