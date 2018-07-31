////////////////////////////////////////////////////////////////////////////////
// Copyright © 2018 Jérémie Dumas
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#include <mma/MMASolver.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>

void ObjSens(double *x, double *x0, int n, int m, double *f, double *df, double *g, double *dg);

// Math helpers
void Print(double *x, int n);
double Min(double d1, double d2);
double Max(double d1, double d2);
int Min(int d1, int d2);
int Max(int d1, int d2);
double Abs(double d1);

int main(int argc, char *argv[]) {
	std::cout<<"///////////////////////////////////////////////////"<<std::endl;
	std::cout<<"Test the MMA algorithm"<<std::endl;
	std::cout<<"///////////////////////////////////////////////////"<<std::endl;

	int dims = 2;
	int npts = 6;
	int n = npts*dims;
	int m = 1;

	double Xmin = 0.0;
	double Xmax = 10.0;
	double movlim = 0.05;

	double *x = new double[n];
	double *xold = new double[n];
	double *x0 = new double[n];
	double *f = new double[1];
	double *df = new double[n];
	double *g = new double[m];
	double *dg = new double[n*m];
	double *xmin = new double[n];
	double *xmax = new double[n];

	std::default_random_engine gen;
	std::uniform_real_distribution<double> dist(Xmin, Xmax);
	for (int i = 0; i < n; ++i) {
		x[i] = dist(gen);
		x0[i] = xold[i] = x[i];
		x[i] -= 1.0;
	}
	Print(x, n);

	for (int i = 0; i < n; ++i) {
		x[i] = dist(gen);
	}
	Print(x, n);

	// Initialize MMA
	MMASolver *mma = new MMASolver(n,m);
	mma->ConstraintModification(true);

	double ch = 1.0;
	int itr = 0;
	while (ch > 0.0002 && itr<1000) {
		itr++;

		ObjSens(x,x0,n,m,f,df,g,dg);

		// Set outer move limits
		for (int i=0;i<n;i++) {
			xmax[i] = Min(Xmax, x[i] + movlim);
			xmin[i] = Max(Xmin, x[i] - movlim);
		}

		// Call the update method
		mma->Update(x,df,g,dg,xmin,xmax);

		// Compute infnorm on design change
		ch = 0.0;
		for (int i=0;i<n;i++) {
			ch = Max(ch,Abs(x[i]-xold[i]));
			xold[i] = x[i];
		}

		// Print to screen
		printf("it.: %d, obj.: %f, ch.: %f \n",itr, f[0], ch);
	}
	Print(x, n);

	// Deallocate
	delete mma;
	return 0;

}

void Print(double *x, int n) {
	std::cout << "x:";
	for (int i=0;i<n;i++) {
		std::cout << " " << x[i];
	}
	std::cout << std::endl;
}

void PrintDg(double *dg, int n, int m, int c) {
	std::cout << "dg_i:";
	for (int i=0;i<n;i++) {
		std::cout << " " << dg[i*m+c];
	}
	std::cout << std::endl;
}

double SquaredNorm(double x, double y) {
	return x*x + y*y;
}

template<typename T, typename U>
double safepow(T x, U n) {
	double y = std::pow(x, n);
	if (std::isinf(y)) {
		return 1e10;
	} else {
		return y;
	}
}

void ObjSens(
	double *x, double *x0, int n, int m, double *f, double *df, double *g, double *dg)
{
	f[0] = 0.0;
	for (int i=0;i<n;i++) {
		f[0] += x[i];
	}

	for (int i=0;i<n;i++) {
		df[i] = 1;
	}

	std::fill_n(dg, n*m, 0.0);

	const int dims = 2;
	const int npts = n / dims;
	const double max_distance = 0.2;
	const double softmax_power = 2.0;

	double s = 0;
	for (int i = 0; i < npts; ++i) {
		double *a = x + i*dims;
		double *b = x0 + i*dims;
		double y[2] = {a[0] - b[0], a[1] - b[1]};
		double r2 = SquaredNorm(y[0], y[1]);
		double dr = softmax_power * safepow(r2, softmax_power - 1.0);
		s += safepow(r2, softmax_power);
		dg[dims*i+0] += dr * 2.0 * y[0];
		dg[dims*i+1] += dr * 2.0 * y[1];
	}
	g[0] = safepow(s / npts, 1.0 / softmax_power) - max_distance;

	for (int i = 0; i < n; ++i) {
		dg[i] *= (1.0 / npts) * (1.0 / softmax_power) * safepow(s / npts, 1.0 / softmax_power - 1.0);
	}

	// for (int k=0;k<q;++k) {
	// 	int i = e[2*k+0];
	// 	int j = e[2*k+1];
	// 	double z[2];
	// 	double r = SpringEnergy(x+2*i, x+2*j, z) - 1.0;
	// 	s += safepow(r, p);
	// 	for (int d = 0; d < 2; ++d) {
	// 		dg[m*(2*i+d)] += p * safepow(r, p - 1.0) * z[d];
	// 		dg[m*(2*j+d)] -= p * safepow(r, p - 1.0) * z[d];
	// 	}
	// }
	// g[0] = safepow(s / q, 1.0 / p) - 0.2;
	// for (int i=0; i<n*m; ++i) {
	// 	dg[i] *= (1.0 / q) * (1.0 / p) * safepow(s / q, 1.0 / p - 1.0);
	// }
	std::cout << g[0] << std::endl;
	//PrintDg(dg, n, m, 0);
}

double Min(double d1, double d2) {
	return d1<d2 ? d1 : d2;
}

double Max(double d1, double d2) {
	return d1>d2 ? d1 : d2;
}

int Min(int d1, int d2) {
	return d1<d2 ? d1 : d2;
}

int Max(int d1, int d2) {
	return d1>d2 ? d1 : d2;
}

double Abs(double d1) {
	return d1>0 ? d1 : -1.0*d1;
}
