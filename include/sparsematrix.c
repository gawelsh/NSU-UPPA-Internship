#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include <assert.h>

#include "sparsematrix.h"
#include "vector.h"


sMatrix createSMatrix(int n, int nnzmax) {
  int i;
  sMatrix m;
  m.n = n;
  m.data = (double *)calloc(nnzmax, sizeof(double));
  m.lstart = (size_t *)calloc(n+1, sizeof(size_t));
  m.j = (size_t *)calloc(nnzmax, sizeof(size_t));
  for (i = 0; i < nnzmax; i++) {
    m.data[i] = 0.0;
    m.j[i] = 0;
  }
  for (i = 0; i < n+1; i++)
    m.lstart[i] = 0.0;
  return m;
}

sMatrix toSMatrix(matrix m) {
  sMatrix sm;
  int i,j, N=0, idx;
  for (i = 0; i < m.n*m.n; i++)
    if(m.data[i]!=0) N++;
  sm = createSMatrix(m.n, N);
  idx = 0;
  for (i = 0; i < m.n; i++) {
    for (j = 0; j < m.n; j++) {
      if(m.data[i*m.n+j]!=0) {
        sm.j[idx] = j;
        sm.data[idx] = m.data[i*m.n+j];
        idx ++;
      }
    }
    sm.lstart[i+1] = idx;
  }
  return sm;
}

sMatrix createSMatrixFromStruct(sMatrix a)  {
  sMatrix b = createSMatrix(a.n, a.lstart[a.n]);
  for(int i=0; i<=a.n; i++)
    b.lstart[i] = a.lstart[i];
  for(int i=0; i<a.lstart[a.n]; i++)
    b.j[i] = a.j[i];
  return b;
}

void destroySMatrix(sMatrix m) {
  m.n = -1;
  free(m.data);
  free(m.j);
  free(m.lstart);
}

void displaySMatrix(sMatrix m) {
  int i,j, nz;
  printf("Sparse matrix %dx%d (%d nz elements):\n",m.n,m.n,(int)m.lstart[m.n]);
  for (i = 0; i < m.n; i++) {
    printf("%3d(%3d)|",i, (int)m.lstart[i]);
    for (nz = 0; nz < m.lstart[i+1]-m.lstart[i]; nz++) {
      j = m.j[m.lstart[i]+nz];
      printf("%8g(%3d) ", m.data[m.lstart[i]+nz], j);
    }
    printf("\n");
  }
}

void prodSMV(sMatrix m, vector a, vector b) {
  assert(a.data!=b.data); // The given vectors must be different
  assert(m.n==a.n && m.n==b.n); // Ensure correctness of elements sizes
  int i,j,k;
  for (i = 0; i < m.n; i++) {
    b.data[i] = 0.0;
    for (k = m.lstart[i]; k < m.lstart[i+1]; k++) {
      j = m.j[k];
      b.data[i] += m.data[k]*a.data[j];
    }
  }
}

void asXpsY(double a, sMatrix X, sMatrix Y, sMatrix res){
	for(int j=0; j<X.lstart[X.n]; j++){
		res.data[j]=a*X.data[j]+Y.data[j];
	}
}

void resetSMatrix(sMatrix m, double a) {
  for(int i=0; i<m.lstart[m.n]; i++)
    m.data[i] = a;
}

void getDiagonal(sMatrix m, vector d) {
  assert(m.n==d.n); // Ensure correctness of elements sizes
  for (int i = 0; i < m.n; i++) {
    for (int nz = 0; nz < m.lstart[i+1]-m.lstart[i]; nz++) {
      int j = m.j[m.lstart[i]+nz];
      if(i==j)
        d.data[i] = m.data[m.lstart[i]+nz];
    }
  }
}
