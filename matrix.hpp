#include <cmath>
#include <climits>
#include <cstdlib>
#include <cstdio>

double*     allocateMatrix      (uint64_t n,uint64_t m) ;
double*     allocateVector      (uint64_t n) ;
void        freeMatrix          (double *A);
void        freeVector          (double *v);
void        setMatrixZero       (double *A, uint64_t n, uint64_t m);
void        setMatrixIdentity   (double *A, uint64_t n);
void        copyMatrix          (double *B, double *A, uint64_t n, uint64_t m) ;
void        writeMatrix         (FILE *stream, double *A, uint64_t n, uint64_t m);
void        absMatrix           (double *Aabs,double *A, uint64_t n, uint64_t m);
double      getMaxInMatrix      (double max, double *A, uint64_t n, uint64_t m);
void        matrixSub           (double *S, double *A, double *B, uint64_t n, uint64_t m);
void        matrixAdd           (double *S, const double *A, const double *B, uint64_t n, uint64_t m);


/**
 * @role create a matrix with random values ranging from -127 to +127, of size size x size.
 * @param size : the lenght of the matrix.
 * @return matrix (size x size).
 */
double * matrixGenerate(uint64_t size);

/**
 * @role : displays a matrix in the form of an array n x n.
 * @param M : the matrix (table).
 * @param n : the length of the matrix (n x n).
 */
void matrixAff(double *M, uint64_t n);

/**
 * @role Performs naive multiplication of matrix A (size p x k) by a matrix B (size k x r).
         The result matrix S = A*B  is of size (k x r).
         We assume that S has already been allocated outside the function.
 * @param S : the final matrix, S = A * B.
 * @param A : the matrix of size p x k.
 * @param B : the matrix of size k x r.
 * @param p : the length of the A matrix.
 * @param k : the width of the A matrix and the length of the B matrix.
 * @param r : the width of the B matrix.
 */
void matrixMultiplyNaive (double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r);

/**
 * @role Performs a multiplication of two square matrices A and B (size n x n) by Strassen algorithm.
         We assume that S has already been allocated outside the function.
 * @param S : the final matrix, S = A * B.
 * @param A : the matrix of size n x n.
 * @param B : the matrix of size n x n.
 * @param n : the length of the matrix.
 */
void matrixMultiplyStrassen (double *S, double *A, double *B, uint64_t n);

/**
 * @role Solves a system of linear equations Ax=b for a double-precision matrix A (size n x n).
         Uses iterative ascension algorithm.
         After the procedure, x contains the solution of Ax=b.
         We assume that x has been allocated outside the function.
 * @param x : the solution Ax=b.
 * @param A : the matrix of size n x n.
 * @param b : the double value.
 * @param n : the length of the matrix.
 */
void SolveTriangularSystemUP (double *x, double *A, double *b, uint64_t n);

/* 
    Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/
bool Triangularize (double *A, double *b, uint64_t n);

/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
        Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool SolveSystemGauss (double *x, double *A, double *b, uint64_t n);
bool TriangularizeLU (double *A, double *L, double *U, uint64_t n);
