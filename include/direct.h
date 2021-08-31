#ifndef DIRECT_H
#define DIRECT_H

#include "vector.h"
#include "matrix.h"


matrix computeLU(matrix A);
void solveLU(matrix LU, vector b, vector x);


#endif /* DIRECT_H */
