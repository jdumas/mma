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

void ObjSens(double *af, double *bf, double *x, int n, int m, double *f, double *df, double *g, double *dg);

// Math helpers
double Min(double d1, double d2);
double Max(double d1, double d2);
int Min(int d1, int d2);
int Max(int d1, int d2);
double Abs(double d1);

int main(int argc, char *argv[]) {
	std::cout<<"///////////////////////////////////////////////////"<<std::endl;
	std::cout<<"Test the MMA algorithm"<<std::endl;
	std::cout<<"///////////////////////////////////////////////////"<<std::endl;

	int n = 100;
	int m = 2;

	double Xmin = -1.0;
	double Xmax = 1.0;
	double movlim = 0.2;

	double *x = new double[n];
	double *xold = new double[n];
	double *f = new double[1];
	double *df = new double[n];
	double *g = new double[m];
	double *dg = new double[n*m];
	double *xmin = new double[n];
	double *xmax = new double[n];

	for (int i=0;i<n;i++) {
		x[i] = 0.5;
		xold[i] = 0.5;
	}

	double *af = new double[n];
	double *bf = new double[n];
	for (int i=0;i<n;i++) {
		af[i] = sin((double)i/((double)(n-1)));
		bf[i] = cos((double)i/((double)(n-1)));
	}

	// Initialize MMA
	MMASolver *mma = new MMASolver(n,m);

	double ch = 1.0;
	int itr = 0;
	while (ch > 0.002 && itr<100) {
		itr++;

		ObjSens(af,bf,x,n,m,f,df,g,dg);

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

	std::cout << "x:";
	for (int i=0;i<n;i++) {
		std::cout << " " << x[i];
	}
	std::cout << std::endl;

	// Deallocate
	delete mma;
	return 0;

}

void ObjSens(
	double *af, double *bf, double *x, int n, int m, double *f, double *df, double *g, double *dg)
{
	f[0] = 0.0;
	for (int i=0;i<n;i++) {
		f[0] = f[0] + pow(x[i],2.0)+ pow(x[i],3.0);
	}

	for (int i=0;i<n;i++) {
		df[i] = 2*x[i]+ 3*pow(x[i],2.0);
		df[i] = df[i]/((double)n);
	}

	g[0] = 0.0;
	for (int i=0;i<n;i++) {
		g[0] += x[i];
	}
	g[0] = g[0]/((double)n) + 1.0;

	g[1] = 0.0;
	for (int i=0;i<n;i++) {
		g[1] += sin(x[i]);
	}
	g[1] = g[1]/((double)n) - 1.0;

	for (int i=0;i<n;i++) {
		dg[i*m+0] = 1.0/((double)n);
		dg[i*m+1] = cos(x[i])/((double)n);
	}
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

