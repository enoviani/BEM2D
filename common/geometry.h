#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

int mod (int a, int b);
void computeNormals(double *x, double *y, int N,double **nx, double **ny, double **H);
void computeMidPoints(double *x, double *y, int N,double **xn, double **yn);
