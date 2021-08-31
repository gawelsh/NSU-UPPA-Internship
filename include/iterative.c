#include "iterative.h"
#include "vector.h"
#include "matrix.h"
#include "sparsematrix.h"
#include <stdio.h>

int solveGaussSeidel(matrix A, vector b, vector x, double tol){
  /* Solves Ax=b with Gauss-Seidel iterative method whithin the given tolerance (tol)
   * Vector x contains x_0 as the initial guess
  * Returns the iteration number.
  */
  
	int k=0;
	int N=1e6;
	int n=b.n;
	double sum;
	vector tampon = createVector(n);
	vector xold = createVectorFromData(x.data,n);
	vector r = createVector(n);

	prodMV(A, x, tampon);
	axpy(-1, tampon, b, r);

	while((normeoo(r)>=tol*normeoo(b))&&(k<N)){
		for(int i=0; i<n; i++){
			sum=0;
			for(int j=0; j<i; j++){
				sum=sum+A.data[i*n+j]*x.data[j];
			}
			for(int j=i+1; j<n; j++){
				sum=sum+A.data[i*n+j]*xold.data[j];
			}
			x.data[i]=(b.data[i]-sum)/A.data[i*n+i];
		}
		prodMV(A, x, tampon);
		axpy(-1, tampon, b, r);
		k=k+1;
		copyVector(xold ,x , n);
	}

	destroyVector(tampon);
	destroyVector(r);
	destroyVector(xold);

	return k;
}

int solveSparseGaussSeidel(sMatrix A, vector b, vector x, double tol){
  /* Solves Ax=b with Gauss-Seidel iterative method whithin the given tolerance (tol)
   * Vector x contains x_0 as the initial guess
  * Returns the iteration number.
  */
  
	int k=0;
	int N=1e6;
	int n=b.n;
	double sum;
	vector diag = createVector(n);
	vector tampon = createVector(n);
	vector xold = createVectorFromData(x.data,n);
	vector r = createVector(n);
	getDiagonal(A, diag);

	prodSMV(A, x, tampon);
	axpy(-1, tampon, b, r);

	while((normeoo(r)>=tol*normeoo(b))&&(k<N)){
		for(int i=0; i<n; i++){
			sum=0;
			for(int j=A.lstart[i]; j<A.lstart[i+1]; j++){
				if(A.j[j]<i){sum+=A.data[j]*x.data[A.j[j]];}
				if(A.j[j]>i){sum+=A.data[j]*xold.data[A.j[j]];}
			}
			x.data[i]=(b.data[i]-sum)/diag.data[i];
		}
		prodSMV(A, x, tampon);
		axpy(-1, tampon, b, r);
		k=k+1;
		copyVector(xold ,x , n);
	}

	destroyVector(tampon);
	destroyVector(r);
	destroyVector(xold);
	destroyVector(diag);

	return k;
}
