#ifndef COMMON_H
#define COMMON_H

#include "matrix.h"

typedef struct {
	double g;     /* Acceleration of gravity */
	double wL;    /* Full wave Length (in m) */
	double genwaves;  /*Number of generated waves*/
	double L;			/* Length (in m) */
	double pi;		/* pi */
	double Hw;		/* Max height of the generated wave (in m)*/
} parameters;

#define PARAMETERS {/*g=*/9.81, /*wL*/2e5,/*genwaves=*/24, /*L=*/5e6, /*pi=*/3.14,/*Hw=*/1}

double Depth(double x);

double N0(double t, parameters p);
double u0(double t, parameters p);

vector initN(vector x, parameters p);
vector initu(vector x, parameters p);
vector initD(vector x, double dx);

#endif /* DIRECT_H */
