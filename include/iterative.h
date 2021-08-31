#ifndef ITERATIVE_H
#define ITERATIVE_H

#include "vector.h"
#include "matrix.h"
#include "sparsematrix.h"

int solveGaussSeidel(matrix A, vector b, vector x, double tol);
int solveSparseGaussSeidel(sMatrix A, vector b, vector x, double tol);

#endif /* ITERATIVE_H */
