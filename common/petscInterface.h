#include<stdio.h>
#include<stdlib.h>
#include <petscksp.h>

void petscSolve(double* matrix,int N, double* rhs, double* x);
void petscMatVecMult(double* matrix,int M, int N, double* vectorx, double* vectory);
void petscInit();
void petscEnd();
