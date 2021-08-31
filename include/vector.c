#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>

#include "vector.h"

vector createVector(int n) {
  /* Create a vector of size n (each element is zero, with calloc) */
  vector u;
  u.n = n;
  // Calloc function is allocating an zero-initialize memory
  u.data = (double *)calloc(n, sizeof(double));
  return u;
}

vector createVectorFromData(const double* const data, int n) {
  /* Creates a vector of size n, with initial data from the given array  */
  vector m = createVector(n);
  for (int i = 0; i < n; i++) {
    m.data[i] = data[i];
  }
  return m;
}

vector linspace(double x0, double xN, int n) {
	/* Create and initialize a vector including the n+1 absisses xi=i*h with h=1/n*/
	vector x=createVector(n+1);
	double h=(xN-x0)/n;
	for(int i=0; i<x.n; i++){
		x.data[i]=x0+i*h;
	}
	return x;
}

void copyVector(vector copy, vector u, int n) {
  /* Copy a vector of size n */
  for (int i = 0; i < n; i++) {
    copy.data[i] = u.data[i];
  }
}

void destroyVector(vector u) {
  /* Deallocate the given vector */
  u.n = -1;
  free(u.data);
}


void displayVector(vector u) {
  /* Display on terminal the given vector */
  int i;
  for (i = 0; i < u.n; i++) {
    printf("%12.6g\n", u.data[i]);
  }
}

double dot(vector a, vector b) {
  double sum=0;

	for(int i=0; i<a.n;i++){
		sum=sum+a.data[i]*b.data[i];
	}

	return sum;
}

double norme1(vector u){
  /* Computes the norm-1 of u :  ||u||_1 */
	double sum=0;
  for(int i=0;i<u.n;i++){
		sum=sum+fabs(u.data[i]);
	}
	return sum;
}

double norme2(vector u){
  /* Computes the norm-2 of u :  ||u||_2 */
	return sqrt(dot(u,u));
}

double normeoo(vector u){
  /* Computes the infinite norm of u : ||u||_oo */
  double c=0;
  for(int i=0;i<u.n;i++){
		if (c<fabs(u.data[i])){
			c=fabs(u.data[i]);
		}
	}
	return c;
}

void axpy(double a, vector X, vector Y, vector res) {
  /* Computes res=a*X+Y with X and Y two given vectors. */
  for(int i=0;i<X.n;i++){
		res.data[i]=a*X.data[i]+Y.data[i];
	}
}

void writeVectorBase(vector **u, int nvec, char fname[50]) {
  /* Write nvec vectors in a file 'fname' vectors elements are placed on the same line. */
  int i,k;
  FILE *f = fopen(fname, "w");
  if (f == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }
  // compute lines number (min of vector lengths)
  int npts = u[0]->n;
  for (k = 1; k < nvec; k++) npts = (npts<u[k]->n)?npts:u[k]->n;
  // Write points
  for (i = 0; i < npts; i++) {
    for (k = 0; k < nvec; k++) {
    fprintf(f, "%g ", u[k]->data[i]);
    }
    fprintf(f, "\n");
  }
  fclose(f);
}

void writeVector(vector u, char fname[50]) {
  vector *uvec = &u;
  writeVectorBase(&uvec, 1, fname);
}
void write2Vector(vector u, vector v, char fname[50]) {
  vector* uvec[2] = {&u, &v};
  writeVectorBase(uvec, 2, fname);
}
void write3Vector(vector u, vector v, vector w, char fname[50]) {
  vector* uvec[3] = {&u, &v, &w};
  writeVectorBase(uvec, 3, fname);
}
void write4Vector(vector u, vector v, vector w, vector y,char fname[50]) {
  vector* uvec[4] = {&u, &v, &w, &y};
  writeVectorBase(uvec, 4, fname);
}
