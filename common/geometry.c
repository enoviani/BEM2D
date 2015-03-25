#include"geometry.h"


int mod (int a, int b)
{
	int ret = a % b;
	if(ret < 0)
		ret+=b;
	return ret;
}

void computeNormals(double *x, double *y, int N,double **nx, double **ny, double **H){
	*nx=malloc(N*sizeof(double));
	*ny=malloc(N*sizeof(double));
	*H=malloc(N*sizeof(double));
	int i,ip;

	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		(*nx)[i]= -(y[ip]-y[i]);
		(*ny)[i]= x[ip]-x[i];
		(*H)[i]= sqrt(pow((*nx)[i],2)+pow((*ny)[i],2));
		(*nx)[i]=(*nx)[i]/(*H)[i];
		(*ny)[i]=(*ny)[i]/(*H)[i];
	}
}

void computeMidPoints(double *x, double *y, int N,double **xn, double **yn){
	*xn=malloc(N*sizeof(double));
	*yn=malloc(N*sizeof(double));
	int i,ip;

	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		(*xn)[i]= (x[ip]+x[i])/2;
		(*yn)[i]= (y[ip]+y[i])/2;
	}
}
