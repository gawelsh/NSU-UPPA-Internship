#include<stdio.h>
#include<math.h>

#include "common.h"
#include "vector.h"

/* Parameters of the wave */

double Depth(double x){
	return 500;
}

double N0(double t, parameters p){
	double tp = p.wL/sqrt(p.g*Depth(0));
	if (t>tp){
		return 0;
	}
	return (p.Hw*(1-cos(2*p.pi*t/tp))/2);
}

double u0(double t, parameters p){
	return N0(t, p)*sqrt(p.g/Depth(0));
}



/* Create first boundary conditions */

vector initN(vector x, parameters p){
	vector N = createVector(x.n);
	N.data[0]=N0(0, p);
	return N;
}

vector initu(vector x, parameters p){
	vector u = createVector(x.n);
	u.data[0]=u0(0, p);
	return u;
}

vector initD(vector x, double dx){
	vector D = createVector(x.n);
	for(int i=0; i<x.n; i++){
		D.data[i]=Depth(i*dx);
	}
	return D;
}
