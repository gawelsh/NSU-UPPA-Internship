#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include "vector.h"
#include "common.h"

parameters p = PARAMETERS;

void timeup(vector N, vector u, vector D, vector x, double t, double dt, double dx) {
	vector stamp_N = createVectorFromData(N.data, N.n);
	vector stamp_u = createVectorFromData(u.data, u.n);

	/* iteration */
	for(int i=1; i<N.n-1; i++){
		N.data[i] -= dt/(2*dx) * ( (D.data[i]+stamp_N.data[i]) * (stamp_u.data[i+1]-stamp_u.data[i-1]) + stamp_u.data[i] * (D.data[i+1]+stamp_N.data[i+1]-D.data[i-1]-stamp_N.data[i-1]) );

		u.data[i] -= dt/(2*dx) * (stamp_u.data[i]*(stamp_u.data[i+1]-stamp_u.data[i-1]) + p.g*(stamp_N.data[i+1]-stamp_N.data[i-1]) );
	}

	/* Dirichlet Boundary condition*/
	N.data[0] = N0(t, p);
	u.data[0] = u0(t, p);

	int i = N.n-1;

// u"(x)=0
	/*N.data[i] -= dt/(2*dx) * ( (D.data[i]+stamp_N.data[i]) * 2*(stamp_u.data[i]-stamp_u.data[i-1]) + stamp_u.data[i] * (D.data[i]+2*stamp_N.data[i]-D.data[i-1]-2*stamp_N.data[i-1]) );

	u.data[i] -= dt/(2*dx) * (stamp_u.data[i]*2*(stamp_u.data[i]-stamp_u.data[i-1]) + p.g*2*(stamp_N.data[i]-stamp_N.data[i-1]) );*/

// Sommerfeld
	double c = sqrt(p.g*D.data[D.n-1]);
	u.data[i] = stamp_u.data[i]-(sqrt(dx*dx+(stamp_u.data[i]-stamp_u.data[i-1])*(stamp_u.data[i]-stamp_u.data[i-1]))/dx)*c*dt/dx*(stamp_u.data[i]-stamp_u.data[i-1]);
	N.data[i] = stamp_N.data[i]-(sqrt(dx*dx+(stamp_u.data[i]-stamp_u.data[i-1])*(stamp_N.data[i]-stamp_N.data[i-1]))/dx)*c*dt/dx*(stamp_N.data[i]-stamp_N.data[i-1]);

	destroyVector(stamp_N);
	destroyVector(stamp_u);
}


int main(int argc, char *argv[argc]){
	printf("Simulation of a tsunami wave. \nNonlinear equation, Euler.\n");

	int n=1000;
	vector x = linspace(0.0, p.L, n);
	double dx = p.L/(n+1);
	double Tps = p.genwaves*p.wL/sqrt(p.g*Depth(0));
	int T=1000;
	double dt=1.0/T;
	double t=0;

	vector N = initN(x, p);
	vector u = initu(x, p);
	vector D = initD(x, dx);

	for(int i=0; i<Tps*T; i++){
		timeup(N, u, D, x, t, dt, dx);
		t += dt;
	}

	char FILE[50];
 	sprintf(FILE, "lw_endNL.dat");
 	write4Vector(x, N, u, D, FILE);

	destroyVector(x);
	destroyVector(N);
	destroyVector(u);
	destroyVector(D);

	return 0;
}
