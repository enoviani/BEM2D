#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include <petscksp.h>
#include"common.h"


void createCircle(int N, double** X, double** Y, double** xn, double** yn, double** nx, double** ny,double** H)
{
	// We first define the points of the boundary

	int i;
	double t;

	*X=malloc(N*sizeof(double));
	*Y=malloc(N*sizeof(double));
	double dt=1.0/N;

	for(i=0;i<N;i++){
		t = 2*M_PI*i*dt;
		(*X)[i]=cos(t);
		(*Y)[i]=sin(t);
	}

	printf("Hi there !\n");
	computeNormals(*X, *Y,N,nx,ny,H);
	printf("Hi there !\n");
	computeMidPoints(*X, *Y, N,xn,yn);
	printf("Hi there !\n");

}

// Validation of the boundary integral problem on the neuman to dirichlet problem on the unit circle
// If we set the neuman data to neum(s)=cos(2*pi*k*s), with s in [0,1[, 
// then we expect to recover the dirichlet data : diri(s)=cos(2*pi*k*s)/abs(k)

void buidInfluenceMatrices(int N, double *X, double *Y, double *xn, double *yn, double *nx, double *ny, double *H, double **M, double** P)
{
	*M=malloc(N*N*sizeof(double));
	*P=malloc(N*N*sizeof(double));

	int i,j,k,jp;

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			k=j+i*N;
			jp=mod(j+1,N);
			if (i==j){
				GreenBoundDiag(H[i],*M+k,*P+k);
			}
			else{
				GreenBound(xn[i],yn[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],*M+k,*P+k);
			}
		}
	}
}

void createPeanut(int N, double** X, double** Y, double** xn, double** yn, double** nx, double** ny,double** H)
{
	// We first define the points of the boundary

	int i;
	double t;

	*X=malloc(N*sizeof(double));
	*Y=malloc(N*sizeof(double));
	double dt=1.0/N;

	for(i=0;i<N;i++){
		t = 2*M_PI*i*dt;
		(*X)[i]=(1-0.5*cos(2*t))*cos(t);
		(*Y)[i]=(1-0.5*cos(2*t))*sin(t);
	}
	computeNormals(*X, *Y,N,nx,ny,H);
	computeMidPoints(*X, *Y, N,xn,yn);

}



int main(int argc,char **argv)
{

	double *X,*Y,*xn,*yn,*nx,*ny,*H;
	// double PI=3.14159265;
	int N;
	int i,j,jp,k;

	if(argc==1)
		N=100;
	else
		N=atoi(argv[1]);

	petscInit();

	createCircle(N,&X,&Y,&xn,&yn,&nx,&ny,&H);
//	createPeanut(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

//initial condition
	double Ux=1;
	double Uy=2;

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	buidInfluenceMatrices(N,X,Y,xn,yn,nx,ny,H,&M,&P);

	// Prepare the Dirichlet data

	double *diri,*neum,*Pneum;
	diri=malloc(N*sizeof(double));
	Pneum=malloc(N*sizeof(double));
	neum=malloc(N*sizeof(double));

	//double wavenumber=5;

	for(i=0;i<N;i++){
	//	neum[i]=cos(2*M_PI*wavenumber*i/N);
		neum[i]=-(nx[i]*Ux+ny[i]*Uy);
	}

	// Here we compute diri=-inv(M)*P*neum
        // first we compute -P*neum
	
	for(i=0;i<N;i++){
		Pneum[i]=-neum[i];
	}
	petscMatVecMult(P,N,N,Pneum,Pneum);
	
	// then, we solve M*diri=P*neum
	
	petscSolve(M, N, Pneum, diri);

	/*double err=0;
	for(i=0;i<N;i++){
		err=err+pow(abs(wavenumber)*diri[i]-neum[i],2);
	}
	err=err/N;
	err=sqrt(err);

	printf("N:\t%i\terr:\t%e\n",N,err);
	*/

//************Building Discretization of domain**********

	int Hx=50;
	int Hy=50;
	int Np=Hx*Hy;

	double mesh1x=-2;
	double mesh2x=2;
	double mesh1y=-2;
	double mesh2y=2;

	double dhx=(mesh2x-mesh1x)/(double)(Hx-1); //Domain x is in [-3,3] with distance between successor and predesessor at x is dhx
	double dhy=(mesh2y-mesh1y)/(double)(Hy-1); //Domain y is in [-3,3] with distance between successor and predesessor at y is dhy

	double* Xgf;
	double* Ygf;
	// Building mesh in array (double*) //
	Xgf=malloc((Np)*sizeof(double));
	Ygf=malloc((Np)*sizeof(double));

	k=0;
	for (i=0;i<Hx;i++){
		k=i*Hx;
		for(j=0;j<Hy;j++){
			Xgf[k+j]=(mesh1x)+dhx*i;
			Ygf[k+j]=(mesh1y)+dhy*j;
		}
	}
//***********************

	// Here, we build the influence matrices M and P for points in (all) domain
	double* Md;
	double* Pd;
	
	double* MDx;
	double* PDx;

	double* MDy;
	double* PDy;


	Md=malloc(Np*N*sizeof(double));
	Pd=malloc(Np*N*sizeof(double));

	MDx=malloc(Np*N*sizeof(double));
	PDx=malloc(Np*N*sizeof(double));

	MDy=malloc(Np*N*sizeof(double));
	PDy=malloc(Np*N*sizeof(double));

	for(i=0;i<Np;i++){
		for(j=0;j<N;j++){
			k=j+i*N;
			jp=mod(j+1,N);
			GreenBound(Xgf[i],Ygf[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],Md+k,Pd+k);
			GreenBoundDx(Xgf[i],Ygf[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],MDx+k,PDx+k);
			GreenBoundDy(Xgf[i],Ygf[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],MDy+k,PDy+k);
		}
	}


	double* Mddiri;
        double* Pdneum;

        double* MDxdiri;
        double* PDxneum;

        double* MDydiri;
        double* PDyneum;

	double* phi;
	double* Vx;
	double* Vy;

        Mddiri=malloc(Np*sizeof(double));
        Pdneum=malloc(Np*sizeof(double));

        MDxdiri=malloc(Np*sizeof(double));
        PDxneum=malloc(Np*sizeof(double));

        MDydiri=malloc(Np*sizeof(double));
        PDyneum=malloc(Np*sizeof(double));

	phi=malloc(Np*sizeof(double));
	Vx=malloc(Np*sizeof(double));
	Vy=malloc(Np*sizeof(double));

	petscMatVecMult(Md, Np, N, diri, Mddiri);
	petscMatVecMult(Pd, Np, N, neum, Pdneum);

	petscMatVecMult(MDx, Np, N, diri, MDxdiri);
	petscMatVecMult(PDx, Np, N, neum, PDxneum);
	
	petscMatVecMult(MDy, Np, N, diri, MDydiri);
	petscMatVecMult(PDy, Np, N, neum, PDyneum);

for(i=0;i<Np;i++){
	phi[i]=Mddiri[i]+Pdneum[i];
        Vx[i]=-MDxdiri[i]+PDxneum[i]-Ux;
        Vy[i]=-MDydiri[i]+PDyneum[i]-Uy;
}

//	displayDomain5(XXV,YYV,VVx,VVy, PPhiv, (Hx*Hy));
	//**** write the result to file domain.txt****
	FILE *fp;

        fp = fopen("domain.txt", "w");

	for(i=0;i<Np;i++){

                fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",Xgf[i],Ygf[i],phi[i],Vx[i],Vy[i]);
        }

	fclose(fp);

	return 0;
}
