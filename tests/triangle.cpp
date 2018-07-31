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

void ObjSens(int *e, int q, double *x, int n, int m, double *f, double *df, double *g, double *dg);

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

	int n = 6*2;
	int m = 1;

	double Xmin = 0.0;
	double Xmax = 10.0;
	double movlim = 0.05;

	double *x = new double[n];
	double *xold = new double[n];
	double *f = new double[1];
	double *df = new double[n];
	double *g = new double[m];
	double *dg = new double[n*m];
	double *xmin = new double[n];
	double *xmax = new double[n];


	// x[2*0+0] = 0; x[2*0+1] = 0;
	// x[2*1+0] = 2; x[2*1+1] = 0;
	// x[2*2+0] = 1; x[2*2+1] = 1;
	// x[2*3+0] = 3; x[2*3+1] = 1;
	// x[2*4+0] = 4; x[2*4+1] = 0;
	// x[2*5+0] = 5; x[2*5+1] = 1;
	{
		int i = 0;
		x[i++] = 0.5;
		x[i++] = 0.0569556;
		x[i++] = 1.5;
		x[i++] = 0.0569556;
		x[i++] = 1;
		x[i++] = 0.922981;
		x[i++] = 3.5;
		x[i++] = 0.967578;
		x[i++] = 4;
		x[i++] = 0.101553;
		x[i++] = 4.5;
		x[i++] = 0.967578;
	}
	for (int i=0;i<n;i++) {
		xold[i] = x[i];
	}
	Print(x, n);

	int q = 6;
	int *e = new int[2*q];
	e[2*0+0] = 0; e[2*0+1] = 1;
	e[2*1+0] = 0; e[2*1+1] = 2;
	e[2*2+0] = 1; e[2*2+1] = 2;
	e[2*3+0] = 3; e[2*3+1] = 4;
	e[2*4+0] = 3; e[2*4+1] = 5;
	e[2*5+0] = 4; e[2*5+1] = 5;

	// Initialize MMA
	MMASolver *mma = new MMASolver(n,m);

	double ch = 1.0;
	int itr = 0;
	while (ch > 0.002 && itr<100) {
		itr++;

		ObjSens(e,q,x,n,m,f,df,g,dg);

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

double Norm(double x, double y) {
	return std::sqrt(x*x + y*y);
}

double SpringEnergy(double *x, double *y, double *z) {
	double r = Norm(x[0]-y[0], x[1]-y[1]);
	z[0] = (x[0]-y[0]) / std::max(r, 1e-6);
	z[1] = (x[1]-y[1]) / std::max(r, 1e-6);
	return r;
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
	int *e, int q, double *x, int n, int m, double *f, double *df, double *g, double *dg)
{
	f[0] = 0.0;
	for (int i=0;i<n;i++) {
		f[0] += x[i];
	}

	for (int i=0;i<n;i++) {
		df[i] = 1;
	}

	std::fill_n(dg, n*m, 0.0);

	const double p = 2;
	double s = 0;
	for (int k=0;k<q;++k) {
		int i = e[2*k+0];
		int j = e[2*k+1];
		double z[2];
		double r = SpringEnergy(x+2*i, x+2*j, z) - 1.0;
		s += safepow(r, p);
		for (int d = 0; d < 2; ++d) {
			dg[m*(2*i+d)] += p * safepow(r, p - 1.0) * z[d];
			dg[m*(2*j+d)] -= p * safepow(r, p - 1.0) * z[d];
		}
	}
	g[0] = safepow(s / q, 1.0 / p) - 0.2;
	for (int i=0; i<n*m; ++i) {
		dg[i] *= (1.0 / q) * (1.0 / p) * safepow(s / q, 1.0 / p - 1.0);
	}
	std::cout << g[0] << std::endl;
	//std::cout << g[1] << std::endl;
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
