#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

void GreenBoundDiag(double H, double* pm, double* pp);
void GreenBound(double xx,double yy,double xm,double ym, double xp,double yp,double nx,double ny,double* pm, double* pp);
void GreenBoundDx(double xx,double yy,double xm,double ym, double xp,double yp,double nx,double ny,double* pm, double* pp);
void GreenBoundDy(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pm, double* pp);

