#include <iostream>
#include "matrix.hpp"

using namespace std;

int main(int argc, char** argv){

    double U[] = {};
    double L[] = {};
    double A[] = {2,2,8,4,3,6,2,4,2};

    TriangularizeLU(A,L,U,3);
    matrixAff(U,3);
    matrixAff(L,3);


    return 0;
}