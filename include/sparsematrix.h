#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <stdlib.h>

#include "vector.h"
#include "matrix.h"

/* Sparse matrix */
typedef struct
{
  int n;        // Matrix size
  double *data; // Matrix data a_{i,j} (size lstart[n])
  size_t *lstart;    // Number of non-zero elements for each line (size n+1)
  size_t *j;       // column index of each non-zero element (size lstart[n])
} sMatrix;


sMatrix createSMatrix(int n, int nnzmax) ;
sMatrix createSMatrixFromStruct(sMatrix a) ;
sMatrix toSMatrix(matrix m);
void destroySMatrix(sMatrix m) ;
void displaySMatrix(sMatrix m) ;
void prodSMV(sMatrix m, vector a, vector b) ;
void asXpsY(double a, sMatrix X, sMatrix Y, sMatrix res) ;
void resetSMatrix(sMatrix m, double a) ;
void getDiagonal(sMatrix m, vector d) ;

#endif /* SPARSEMATRIX_H */
