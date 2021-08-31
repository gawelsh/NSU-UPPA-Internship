#ifndef VECTOR_H
#define VECTOR_H

/* Vector */
typedef struct
{
  int n;        // Vector length
  double *data; // Vector data
} vector;

// Create and destroy vectors
vector createVector(int n) ;
vector createVectorFromData(const double * const data, int n) ;
vector linspace(double x0, double xN, int n) ;
void copyVector(vector copy, vector u, int n) ;
void destroyVector(vector v) ;

// Display vectors (only for small vectors)
void displayVector(vector m) ;

// compute vector norms
double norme1(vector u) ;
double norme2(vector u) ;
double normeoo(vector u) ;

// Compute ax+y
void axpy(double a, vector x, vector y, vector res);
// compute dot product
double dot(vector a, vector b);

// Save vector values in files
void writeVectorBase(vector **u, int nvec, char fname[50]);
void writeVector(vector u, char fname[50]) ;
void write2Vector(vector u, vector v, char fname[50]) ;
void write3Vector(vector u, vector v, vector w, char fname[50]) ;
void write4Vector(vector u, vector v, vector w, vector y, char fname[50]) ;

#endif /* VECTOR_H */
