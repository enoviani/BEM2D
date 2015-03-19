#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"green.h"
#include<complex.h>
#include"petscInterface.h"
#include <petscksp.h>


int mod (int a, int b)
{
	int ret = a % b;
	if(ret < 0)
		ret+=b;
	return ret;
}

void createCircle(int N, double** X, double** Y, double** xn, double** yn, double** nx, double** ny,double** H)
{

	// We first define the points of the boundary

	double PI=3.14159265;
	int i,ip;
	double t;

	*X=malloc(N*sizeof(double));
	*Y=malloc(N*sizeof(double));
	double dt=1.0/N;

	for(i=0;i<N;i++){
		t = 2*PI*i*dt;
		(*X)[i]=cos(t);
		(*Y)[i]=sin(t);
	}

	// Then we compute the segment middle points

	*xn = malloc(N*sizeof(double));
	*yn = malloc(N*sizeof(double));

	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		(*xn)[i]= ((*X)[ip]+(*X)[i])/2;
		(*yn)[i]= ((*Y)[ip]+(*Y)[i])/2;
	}

	// Normal vectors
	
	*nx=malloc(N*sizeof(double));
	*ny=malloc(N*sizeof(double));
	*H=malloc(N*sizeof(double));

	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		(*nx)[i]= -((*Y)[ip]-(*Y)[i]);
		(*ny)[i]= (*X)[ip]-(*X)[i];
		(*H)[i]= sqrt(pow((*nx)[i],2)+pow((*ny)[i],2));
		(*nx)[i]=(*nx)[i]/(*H)[i];
		(*ny)[i]=(*ny)[i]/(*H)[i];
	}
}

// Validation of the boundary integral problem on the neuman to dirichlet problem on the unit circle
// If we set the neuman data to neum(s)=cos(2*pi*k*s), with s in [0,1[, 
// then we expect to recover the dirichlet data : diri(s)=cos(2*pi*k*s)/abs(k)

int main(int argc,char **argv)
{

	double *X,*Y,*xn,*yn,*nx,*ny,*H;
	double PI=3.14159265;
	int N;
	int i,j,jp,k;

	if(argc==1)
		N=10;
	else
		N=atoi(argv[1]);

	petscInit();

	createCircle(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	M=malloc(N*N*sizeof(double));
	P=malloc(N*N*sizeof(double));

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			k=i+j*N;
			jp=mod(j+1,N);
			if (i==j){
				M[k]=0.5;
				P[k]=(1/(2*PI))*((H[i]*log(H[i]/2))-H[i]);
			}
			else{
				GreenBound(xn[i],yn[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],&(M[k]),&(P[k]));
			}
		}
	}

	// Prepare the Dirichlet data

	double *diri,*neum,*Pneum;
	diri=malloc(N*sizeof(double));
	Pneum=malloc(N*sizeof(double));
	neum=malloc(N*sizeof(double));

	double wavenumber=5;

	for(i=0;i<N;i++){
		neum[i]=cos(2*PI*wavenumber*i/N);
	}

	// Here we compute diri=-inv(M)*P*neum
        // first we compute -P*neum
	
	for(i=0;i<N;i++){
		Pneum[i]=-neum[i];
	}
	petscMatVecMult(P,N,N,Pneum,Pneum);
	
	// then, we solve M*diri=P*neum
	
	petscSolve(M, N, Pneum, diri);

	double err=0;
	for(i=0;i<N;i++){
		err=err+pow(abs(wavenumber)*diri[i]-neum[i],2);
//		printf("%f \t %f \n",abs(wavenumber)*diri[i],neum[i]);
	}
	err=err/N;
	err=sqrt(err);

	printf("Number of points: %i \t Error: %e\n",N,err);


	return 0;
}
