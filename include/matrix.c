#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<assert.h>

#include "matrix.h"
#include "vector.h"

matrix createMatrix(int n) {
  /* Creates a square matrix of size nxn */
  int i;
  matrix m;
  m.n = n;
  m.data = (double *)calloc(n*n, sizeof(double));
  for (i = 0; i < m.n*m.n; i++)
    m.data[i] = 0.0;
  return m;
}

matrix createMatrixFromData(double *data, int n) {
  /* Creates a matrix of size nxn, with initial data from the given array  */
  matrix m;
  m.n = n;
  m.data = (double *)calloc(n*n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      m.data[i*n+j] = data[i*n+j];
    }
  }
  return m;
}

void destroyMatrix(matrix m) {
  /* Free the memory for this matrix */
  m.n = -1;
  free(m.data);
}

void displayMatrix(matrix m) {
  /* Display the given matrix (use only for small sizes) */
  int i,j;
  printf("[");
  for (i = 0; i < m.n; i++) {
    printf("%s[",((i==0)?"":" "));
    for (j = 0; j < m.n; j++) {
      printf("%12.6g ", m.data[i*m.n+j]);
      if(j < m.n-1)
        printf(",");
    }
    printf("]");
    if(i < m.n-1)
      printf(",\n");
  }
  printf("]\n");
}

void aXpY(double a, matrix X, matrix Y, matrix res) {
	for(int i=0; i<X.n; i++){
		for(int j=0; j<X.n; j++){
			res.data[i*res.n+j]=a*X.data[i*X.n+j]+Y.data[i*Y.n+j];
		}
	}
}

void prodMV(matrix m, vector a, vector b) {
  /* Computes the b=Ma matrix-vector product */
  int i,j;

	for(i=0;i<m.n;i++){
		b.data[i]=0;
	}

	for(i=0; i<a.n; i++) {
		for(j=0; j<a.n; j++) {
			b.data[i]=b.data[i]+m.data[i*m.n+j]*a.data[j];
		}
	}
}
