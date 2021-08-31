#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include "vector.h"
#include "matrix.h"
#include "iterative.h"
#include "common.h"

parameters p = PARAMETERS;

sMatrix initA(vector u, vector D, double dt, double dx) {
	sMatrix A=createSMatrix(u.n, 14+5*(u.n-4));
	A.lstart[0]=0;
	double K = p.g*dt*dt/16/dx/dx;

	// Left Boundary Conditions
	A.data[0]=1; A.j[0]=0;
	A.data[1]=0; A.j[1]=1;
	A.data[2]=0; A.j[2]=2;
	A.lstart[1]=3;

	A.data[3]=0; A.j[3]=0;
	A.data[4]=1-D.data[2]*K; A.j[4]=1;
	A.data[5]=(D.data[1]-D.data[3])*K; A.j[5]=2;
	A.data[6]=-D.data[2]*K; A.j[6]=3;
	A.lstart[2]=7;

  // Main equation
	int k=7; int l=3;
  for(int i=2; i<u.n-2; i++){
		A.data[k]=-D.data[i-1]*K;		A.j[k]=i-2;
		A.data[k+1]=(D.data[i]-D.data[i-2])*K; 	A.j[k+1]=i-1;
		A.data[k+2]=1+(D.data[i+1]+D.data[i-1])*K; 	A.j[k+2]=i;
		A.data[k+3]=(D.data[i]-D.data[i+2])*K; 	A.j[k+3]=i+1;
		A.data[k+4]=-D.data[i+1]*K; 	A.j[k+4]=i+2;
		k=k+5;
		A.lstart[l]=A.lstart[l-1]+5;	l++;
	}

	// Right Boundary Condition
	A.data[k]=-D.data[D.n-3]*K; 	A.j[k]=u.n-4;
	A.data[k+1]=(D.data[D.n-2]-D.data[D.n-4])*K; 	A.j[k+1]=u.n-3;
	A.data[k+2]=1+D.data[D.n-3]*K; A.j[k+2]=u.n-2;
	A.data[k+3]=0; 	A.j[k+3]=u.n-1;
	A.lstart[l]=A.lstart[l-1]+4; l++;

	A.data[k+4]=0; 	A.j[k+4]=u.n-3;
	A.data[k+5]=0; 	A.j[k+5]=u.n-2;
	A.data[k+6]=1; 	A.j[k+6]=u.n-1;
	A.lstart[l]=A.lstart[l-1]+3; //Sommerfeld 

	/* A.data[k]=-K*Depth(1); 	A.j[k]=u.n-4;
	A.data[k+1]=K*(Depth(1)-Depth(1)); 	A.j[k+1]=u.n-3;
	A.data[k+2]=1+K*(2*Depth(1)+Depth(1));	A.j[k+2]=u.n-2;
	A.data[k+3]=-K*(3*Depth(1)-Depth(1)); 	A.j[k+3]=u.n-1;
	A.lstart[l]=A.lstart[l-1]+4; l++;

	A.data[k+4]=-K*Depth(1); 	A.j[k+4]=u.n-3;
	A.data[k+5]=K*(3*Depth(1)-Depth(1)); 	A.j[k+5]=u.n-2;
	A.data[k+6]=1+K*(3*Depth(1)-Depth(1)); 	A.j[k+6]=u.n-1;
	A.lstart[l]=A.lstart[l-1]+3; //u"(x)=0 */

  return A;
}

vector initb(vector u, vector N, vector D, double t, double dt, double dx) {
  vector b;
  b = createVector(u.n);
	double C = dt/4/dx;

	// Left Boundary Conditions	
	b.data[0]=N0(t, p)*sqrt(p.g/D.data[0]);
	b.data[1]=u.data[1]-p.g*C*(2*N.data[2]-N.data[0]-N0(t+dt, p)-C*((D.data[3]-D.data[1])*u.data[2]+D.data[2]*(u.data[3]-u.data[1])));

	// Main Equation
	for(int i=2; i<u.n-2; i++){
		b.data[i]=u.data[i]-p.g*C*(2*N.data[i+1]-2*N.data[i-1]-C* ( (D.data[i+2]-D.data[i])*u.data[i+1] + D.data[i+1]*(u.data[i+2]-u.data[i]) - (D.data[i]-D.data[i-2])*u.data[i-1] - D.data[i-1]*(u.data[i]-u.data[i-2]) ) ) ;
	}

	// Right Boundary Conditions
//Sommerfeld
	double c=(sqrt(dx*dx+(u.data[u.n-1]-u.data[u.n-2])*(u.data[u.n-1]-u.data[u.n-2]))/dx)*sqrt(D.data[D.n-1]*p.g);

	b.data[u.n-2]=u.data[u.n-2]-p.g*C* (2*N.data[u.n-1]-2*N.data[u.n-3] - c*dt/dx*(N.data[u.n-1]-N.data[u.n-2]) + C*((D.data[D.n-2]-D.data[D.n-4])*u.data[u.n-3] - D.data[D.n-3]*(u.data[u.n-2]-u.data[u.n-4]) ) );

	b.data[u.n-1]=u.data[u.n-1]-c*dt/dx*(u.data[u.n-1]-u.data[u.n-2]);

//u"(L)=0
	/*b.data[u.n-1]=u.data[u.n-1]-2*p.g*C* (2*N.data[u.n-1]-2*N.data[u.n-2]-C* ( (Depth(1)-Depth(1))*u.data[u.n-1] + 2*Depth(1)*(u.data[u.n-1]-u.data[u.n-2]) - (Depth(1)-Depth(1))*u.data[u.n-2] - Depth(1)*(u.data[u.n-1]-u.data[u.n-3]) ) );

	b.data[u.n-2]=u.data[u.n-2]-p.g*C* (2*N.data[u.n-1]-2*N.data[u.n-3]-C* ( (Depth(1)-Depth(1))*u.data[u.n-1] + 2*Depth(1)*(u.data[u.n-1]-u.data[u.n-2]) - (Depth(1)-Depth(1))*u.data[u.n-3] - Depth(1)*(u.data[u.n-2]-u.data[u.n-4]) ) );*/

  return b;
}

void timeup(vector N, vector u, vector D, vector x, double t, double dt, double dx) {
	sMatrix A = initA(u, D, dt, dx);
	vector b = initb(u, N, D, t, dt, dx);
	vector stamp_u = createVectorFromData(u.data, u.n);
	vector stamp_N = createVectorFromData(N.data, N.n);
	
	solveSparseGaussSeidel(A, b, u, 1e-8);


	for(int i=1; i<N.n-1; i++){
		N.data[i] -= dt/4/dx*((D.data[i+1]-D.data[i-1])*(u.data[i]+stamp_u.data[i]) + D.data[i]*(u.data[i+1]-u.data[i-1]+stamp_u.data[i+1]-stamp_u.data[i-1]));
	}

	// Boundary Conditions
	N.data[0]=N0(t, p);

	double c = (sqrt(dx*dx+(stamp_u.data[u.n-1]-stamp_u.data[u.n-2])*(stamp_u.data[u.n-1]-stamp_u.data[u.n-2]))/dx)*sqrt(p.g*D.data[D.n-1]);

	N.data[N.n-1] -= c*dt/dx*(stamp_N.data[N.n-1]-stamp_N.data[N.n-2]);  // Sommerfeld's Boundary Condition

	/* N.data[N.n-1] -= dt/4/dx*((D.data[D.n-1]-D.data[D.n-2])*(u.data[N.n-1]+stamp_u.data[N.n-1]) + 2*D.data[D.n-1]*(u.data[N.n-1]-u.data[N.n-2]+stamp_u.data[N.n-1]-stamp_u.data[N.n-2])); //u"(x)=0 */

	destroySMatrix(A);
	destroyVector(b);
	destroyVector(stamp_u);
	destroyVector(stamp_N);
}


int main(int argc, char *argv[argc]){
	printf("Simulation of a tsunami wave. \nLinear equation, Crank-Nicolson. \n");

	int n=1000;
	vector x = linspace(0.0, p.L, n);
	double dx = p.L/(n+1);
	double Tps = p.genwaves*p.wL/sqrt(p.g*Depth(0));
	int T=10000;
	double dt=1.0/T;
	double t=0;

	vector N = initN(x, p);
	vector u = initu(x, p);
	vector D = initD(x, dx);

	for(int i=0; i<Tps*T; i++){
		timeup(N, u, D, x, t, dt, dx);
		t+=dt;
	}

	char FILE[50];
 	sprintf(FILE, "lw_end.dat");
 	write4Vector(x, N, u, D, FILE);

	destroyVector(x);
	destroyVector(N);
	destroyVector(u);
	destroyVector(D);

	return 0;
}
