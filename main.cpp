/*
 * @nameProject :\n matrix library.
 *
 * @project :\n The aim of this project is to create a matrix library, in order to be able to reuse
 *            it in other projects later.
 *
 * @partOfTheProject :\n project as part of courses at the University of Nantes.
 *
 * @lastUpdate :\n 20 february by Bastien COUTAND and Cyprien GARNIER.
 *
 * @Creator :\n COUTAND Bastien (bastien.coutand@etu.univ-nantes.fr),
 *            GARNIER Cyprien (cyprien.garnier@etu.univ-nantes.fr).
 */
#include <iostream>
#include "matrix.hpp"

using namespace std;

int main(int argc, char** argv) {

    


    double L[] = {
            1, 0,  0,  0,
            6, 1, 0,  0,
            2, 19, 1, 0,
            4, 10, 11, 1
    };


    double U[] = {
           10, 3, 8, 7,
           0,  5, 3, 5,
           0,  0, 2, 3,
           0,  0, 0, 6
    };

    matrixAff(L,4,4);
    matrixAff(U,4,4);

    double *E = fusionLU(L,U,4);

    matrixAff(E,4,4);
    


    return 0;
}