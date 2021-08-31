#include <math.h>
#include <stdlib.h>

#include "vector.h"
#include "matrix.h"
#include "direct.h"

matrix computeLU(matrix A) {
  /* Compute LU decomposition of given matrix A
   * Result L and U are stored as a square matrix (without L diagonal, which is 1)
   *
   * Result : ( U  U  U  U )
   *          ( L  U  U  U )
   *          ( L  L  U  U )
   *          ( L  L  L  U )
   */
  matrix lu;
  int i,j,r;
  // Copy A into LU matrix
  lu = createMatrixFromData(A.data, A.n);
  // Compute LU decomposition
  for (r=0; r<A.n-1; r++) {
    for (i=r+1; i<A.n; i++) {
      lu.data[i*A.n+r] = lu.data[i*A.n+r]/lu.data[r*A.n+r];
      for (j=r+1; j<A.n; j++) {
        lu.data[i*A.n+j] -= lu.data[i*A.n+r]*lu.data[r*A.n+j];
      }
    }
  }
  return lu;
}

void forward_subst(matrix L, vector b, vector y) {
  int i,j;
  int n = L.n;
  for (i = 0; i < n; i++) {
    y.data[i] = b.data[i];
    for (j=0; j<i; j++) {
      y.data[i] -= L.data[i*n+j]*y.data[j];
    }
  }
}

void backward_subst(matrix U, vector y, vector x) {
  int i,j;
  int n = U.n;
  for (i=n-1; i>=0; i--) {
    x.data[i] = y.data[i];
    for (j=i+1; j<n; j++) {
      x.data[i] -= U.data[i*n+j]*x.data[j];
    }
    x.data[i] /= U.data[i*n+i];
  }
}

void solveLU(matrix LU, vector b, vector x){
	vector y = createVector(b.n);

	LU = computeLU(LU);

	forward_subst(LU,b,y);
	backward_subst(LU,y,x);

	destroyVector(y);
}
