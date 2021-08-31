#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

/* Square matrix */
typedef struct
{
  int n;        // Matrix size
  double *data; // Matrix data a_{i,j} = data[i*n+j]
} matrix;


// Create and destroy matrix
matrix createMatrix(int n) ;
matrix createMatrixFromData(double *data, int n) ;
void destroyMatrix(matrix m) ;

// Display matrix (only for small matrices)
void displayMatrix(matrix m) ;

// Add a*X to Y
void aXpY(double a, matrix X, matrix Y, matrix res);

// Matrix-vector product
void prodMV(matrix m, vector a, vector b) ;

#endif /* MATRIX_H */
