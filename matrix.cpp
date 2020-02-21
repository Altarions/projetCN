#include <cstdint>
#include <iostream>
#include <stdlib.h>
#include "matrix.hpp"

using namespace std;

/**
 * @role memory allocation for a matrix of size n x m and initialization to 0.
 * @param n : the length of the matrix.
 * @param m : the width of the matrix.
 * @return matrix (n x m).
 */
double *allocateMatrix(uint64_t n,uint64_t m) {
    double *A;
    A = (double *) calloc (n * m, sizeof(double));
    return A;
}

/**
 * @role frees the memory allocated to matrix A.
 * @param A : the matrix.
 */
void freeMatrix(double *A) {
    free(A);
}

/**
 * @role allocates a n sized vector and initializes all entries to 0.
 * @param n : the length of the matrix.
 * @return vector of size n x n.
 */
double *allocateVector(uint64_t n) {
    double *v;
    v = (double *) calloc(n, sizeof(double));
    return v;
}

/**
 * @role : trees the memory allocated to a vector.
 * @param v : the vector.
 */
void freeVector(double *v) {
    free(v);
}

/**
 * @role sets a n * m matrix A to all zeros.
 * @param A : the matrix that we want to set to 0.
 * @param n : the length of the matrix.
 * @param m : the width of the matrix.
 */
void setMatrixZero(double *A, uint64_t n, uint64_t m) {
    uint64_t i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
        /* Note that for a n x m matrix flattened to a 1D array, 
        element A_ij has index i * m + j
        */
        A[i * m + j] = 0.0;
        }
    }
}

/**
 * @role sets a n * n matrix A to identity.
 * @param A : the identity matrix.
 * @param n : the lenght of the matrix.
 */
void setMatrixIdentity (double *A, uint64_t n) {
    uint64_t i, j;

    for (i = 0; i < n; i++) {
        for (j = 0;j < n; j++) {
            A[i * n + j] = 0.0;
        }
        A[i * n + i] = 1.0;
    }
}

/**
 * @role copies a matrix.
 * @param B : the matrix of size n x n.
 * @param A : the matrix copied of size m x m.
 * @param n : the length of the matrix.
 * @param m : the width of the matrix.
 */
void copyMatrix (double *B, const double *A, uint64_t n, uint64_t m) {
    uint64_t i,j;

    for (i = 0; i < n ; i++) {
        for (j = 0; j < n; j++) {
            B[i * m + j] = A[i * n + j];
        }
    }
}

/**
 * @role Writes a matrix to a stream. For example, writing a matrix to standard output is
         writeMatrix(stdout, A, n, m);
         A sream can also be a file.
 * @param stream : the place where we want to write.
 * @param A : the matrix to write, of size n x m.
 * @param n : the length of the matrix.
 * @param m : the width of the matrix.
 */
void writeMatrix (FILE *stream, double *A, uint64_t n, uint64_t m) {
	fprintf(stream, "%d %d \n", (int)n, (int)m);
	int i, j;

	for (i = 0; i < n; ++i) {
	      for (j = 0; j < m; ++j) {
		      fprintf(stream, "%f \t", A[i * m + j]);
	      }
	      fprintf(stream, "\n");
	}
}

/**
 * @role the function computes the element-by-element abs of matrix A.
 * @param A : the matrix of size n x m.
 * @param n : the length of the matrix.
 * @param m : the width of the matrix.
 */
void absMatrix (double *A, uint64_t n, uint64_t m) {
	uint64_t i,j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
            A[i*m + j] = fabs(A[i*m + j]);
		}
	}
}

/**
 * @role Performs addition of two matrix A (size n x m) and B (size n x m).
         The result S = A + B is a n x m matrix.
         We consider that S is allocated outside the function.
 * @param S : the final matrix, S = A + B.
 * @param A : the matrix of size n x n.
 * @param B : the matrix of size m x m.
 * @param n : the length of the A matrix.
 * @param m : the length of the B matrix.
 */
void matrixAdd(double *S, const double *A, const double *B, uint64_t n, uint64_t m){
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] + B[i*m + j];
		}
	}
}

/**
 * @role Performs subtraction of two matrix A (size n x m) and B (size n x m).
         The result S = A - B is a n x m matrix.
         We consider that S is allocated outside the function.
 * @param S : the final matrix, S = A - B.
 * @param A : the matrix of size n x n.
 * @param B : the matrix of size m x m.
 * @param n : the length of the A matrix.
 * @param m : the length of the B matrix.
 */
void matrixSub(double *S, double *A, double *B, uint64_t n, uint64_t m){
    uint64_t i,j;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < m; ++j)
		{
            S[i*m + j] = A[i*m + j] - B[i*m + j];
		}
	}
}

/**
 * @role For a double m x n matrix A the function returns its maximum in absolute value
         element.
 * @param max : the maximum values in the matrix.
 * @param A : the matrix of size n x m.
 * @param n : the length of the matrix.
 * @param m : the width of the matrix.
 * @return the maximum values in the A matrix.
 */
double getMaxInMatrix(double max, double *A, uint64_t n, uint64_t m) {
	double maxA = fabs(A[0]);
	double current = fabs(A[0]);
	int i,j;

	for(i = 0; i < n; ++i) {
		for(j = 0; j < m; ++j) {
			current = fabs(A[i * m + j]);
			if(current > maxA) maxA = current;
		}
	}
    return maxA;
}

/**
 * @role : displays a matrix in the form of an array n x n.
 * @param M : the matrix (table).
 * @param n : the length of the matrix (n x n).
 */
void matrixAff(double *M, uint64_t n) {
    cout << "\nLa matrice est de taille " << n << " et de forme :" << endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << M[i * n + j] << "\t";
        }
        cout << "\n";
    }
}

/**
 * @role create a matrix with random values ranging from -127 to +127, of size size x size.
 * @param size : the lenght of the matrix.
 * @return matrix (size x size).
 */
double * matrixGenerate(uint64_t size) {
    double *A = allocateMatrix(size, size);
    double randomNumber = 0;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A[i * size + j] = rand() % (127 - (-127) + 1) - 127;
        }
    }
    // matrixAff(A, size); affichage de la nouvelle matrice.
    return A;
}


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
void matrixMultiplyNaive (double *S, double *A, double *B, uint64_t p, uint64_t k, uint64_t r){
    for(int i=0; i<k;i++){
        for(int j=0; j<r; j++){
            S[i*p+j]=0;
            for(int m=0; m<k; m++){
                S[i*p+j] += A[i*k+m] * B[m*r+j];
            }
        }
    }
}

/**
 * @role Performs a multiplication of two square matrices A and B (size n x n) by Strassen algorithm.
         We assume that S has already been allocated outside the function.
 * @param S : the final matrix, S = A * B.
 * @param A : the matrix of size n x n.
 * @param B : the matrix of size n x n.
 * @param n : the length of the matrix.
 */
void matrixMultiplyStrassen (double *S, double *A, double *B, uint64_t n) {

    /* Votre code ici */
}

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
void SolveTriangularSystemUP (double *x, double *A, double *b, uint64_t n) {

    /* Votre code ici */
}

/* 
    Performs Gauss elimination for given a matrix A (size n x n) and a vector b (size n).
    Modifies directly matrix A and vector b.
    In the end of the procedure, A is upper truangular and b is modified accordingly.
    Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is impossible to triangularize. 
*/
bool Triangularize (double *A, double *b, uint64_t n) {
    
    

    for(int cl = 0; cl<n-1; cl++){
        for(int lgn = cl+1; lgn<n; lgn++){
            if(A[lgn*n+cl] != 0){
                int x= A[lgn*n+cl]/A[cl];
                for(int w=0;w<n;w++){
                    A[lgn*n+w]=A[lgn*n+w] - x*A[w];
                }
                b[lgn]-= x*b[lgn];
            }
        }
    }
    for(int i = 0; i < n; ++i){
        if (A[i*n + i] == 0){
            return false;
        }
	}

    matrixAff(A,n);
    return true;
}

/*
    Solves a system of linear equations Ax=b, given a matrix A (size n x n) and vector b(size n).
    Uses Gauss elimination algorithm based on truangularization and the ascension solving.
    After the procedure, vector x contains the solution to Ax=b.
    We assume that x has been allocated outside the function.
        Returns a boolean variable: 
        *  true in case of success and 
        *  false in case of failure, for example matrix is of rank <n .
*/
bool SolveSystemGauss (double *x, double *A, double *b, uint64_t n) {
    
    /* Votre code ici */

    return false;
}

bool TriangularizeLU (double *A, double *L, double *U, uint64_t n) {
    
    U = A;
    setMatrixIdentity(L,n);

    for(int cl = 0; cl<n-1; cl++){
        for(int lgn = cl+1; lgn<n; lgn++){
            if(U[lgn*n+cl] != 0){
                int x= U[lgn*n+cl]/U[cl];
                for(int w=0;w<n;w++){
                    U[lgn*n+w]=U[lgn*n+w] - x*U[w];
                }
                U[lgn*n+cl]=x;
            }
        }
    }
    for(int i = 0; i < n; ++i){
        if (A[i*n + i] == 0){
            return false;
        }
	}

    matrixAff(A,n);
    return true;
}
double * fusionLU (double *L, double *U, uint64_t n){
    double *fus=allocateMatrix(n,n);
    for(int cl = 0; cl<n-1; cl++){
        for(int lgn = cl+1; lgn<n; lgn++){
            
        }
    }
}