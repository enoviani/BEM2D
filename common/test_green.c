#include"common.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


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
		(*X)[i]=(1+0.5*cos(2*t))*cos(t);
		(*Y)[i]=(1+0.5*cos(2*t))*sin(t);
	}
	computeNormals(*X, *Y,N,nx,ny,H);
	computeMidPoints(*X, *Y, N,xn,yn);

}


// The peanut testcase for Dirichlet to Neumann (see the doc)
//
double test_peanutDN(int N)
{
	double *X,*Y,*xn,*yn,*nx,*ny,*H;
	// double PI=3.14159265;
	int i;

	createPeanut(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	buidInfluenceMatrices(N,X,Y,xn,yn,nx,ny,H,&M,&P);

	// Prepare the Dirichlet data

	double *diri,*neum,*rhs,*neumex;
	diri=malloc(N*sizeof(double));
	rhs=malloc(N*sizeof(double));
	neum=malloc(N*sizeof(double));
	neumex=malloc(N*sizeof(double));

	for(i=0;i<N;i++){
		diri[i]=(1/(2*M_PI))*log(1.0/sqrt(xn[i]*xn[i]+yn[i]*yn[i]));
		neumex[i]=-(1/(2*M_PI))*(1.0/(xn[i]*xn[i]+yn[i]*yn[i]))*(xn[i]*nx[i]+yn[i]*ny[i]);
	}

	// Here we compute neum=- P \(M*(diri'));
        // first we compute rhs=-M*diri
	
	for(i=0;i<N;i++){
		rhs[i]=-diri[i];
	}
	petscMatVecMult(M,N,N,rhs,rhs);
	
	// then, we solve P*neum=rhs
	
	petscSolve(P, N, rhs, neum);

	double err=0;
	for(i=0;i<N;i++){
		err=err+pow(neumex[i]-neum[i],2);
	}
	
	err=err/N;
	err=sqrt(err);
	return err;
}

// The peanut testcase for Neumann to Dirichlet (see the doc)
//
double test_peanutND(int N)
{
	double *X,*Y,*xn,*yn,*nx,*ny,*H;
	// double PI=3.14159265;
	int i;

	createPeanut(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	buidInfluenceMatrices(N,X,Y,xn,yn,nx,ny,H,&M,&P);

	// Prepare the Dirichlet data

	double *diri,*neum,*rhs,*diriex;
	diri=malloc(N*sizeof(double));
	rhs=malloc(N*sizeof(double));
	neum=malloc(N*sizeof(double));
	diriex=malloc(N*sizeof(double));

	for(i=0;i<N;i++){
		diriex[i]=(1/(2*M_PI))*log(1.0/sqrt(xn[i]*xn[i]+yn[i]*yn[i]));
		neum[i]=-(1/(2*M_PI))*(1.0/(xn[i]*xn[i]+yn[i]*yn[i]))*(xn[i]*nx[i]+yn[i]*ny[i]);
	}

	// Here we compute diri=- M \(P*(neum));
        // first we compute rhs=-P*neum
	
	for(i=0;i<N;i++){
		rhs[i]=-neum[i];
	}
	petscMatVecMult(P,N,N,rhs,rhs);
	
	// then, we solve M*diri=rhs
	
	petscSolve(M, N, rhs, diri);

	double err=0;
	for(i=0;i<N;i++){
		err=err+pow(diriex[i]-diri[i],2);
	}
	
	err=err/N;
	err=sqrt(err);
	return err;
}

double test_peanutBulk(int N)
{

	// geometric data of the boundary
	double *X,*Y,*xn,*yn,*nx,*ny,*H;


	// points inside the bulk
	
	double Xb[4]={2,2,-2,-2};
	double Yb[4]={2,-2,2,-2};

	int i,j,k,jp;

	createPeanut(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	M=malloc(4*N*sizeof(double));
	P=malloc(4*N*sizeof(double));

	for(i=0;i<4;i++){
		for(j=0;j<N;j++){
			k=j+i*N;
			jp=mod(j+1,N);
			GreenBound(Xb[i],Yb[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],M+k,P+k);
		}
	}

	// Prepare the Dirichlet and Neumann data, compute the exact solution.

	double *neumex,*diriex;
	double *Pneumex,*Mdiriex;
	double *solex,*sol;

	Mdiriex=malloc(4*sizeof(double));
	diriex=malloc(N*sizeof(double));
	Pneumex=malloc(4*sizeof(double));
	neumex=malloc(N*sizeof(double));
	solex=malloc(4*sizeof(double));
	sol=malloc(4*sizeof(double));

	for(i=0;i<N;i++){
		diriex[i]=(1/(2*M_PI))*log(1.0/sqrt(xn[i]*xn[i]+yn[i]*yn[i]));
		neumex[i]=-(1/(2*M_PI))*(1.0/(xn[i]*xn[i]+yn[i]*yn[i]))*(xn[i]*nx[i]+yn[i]*ny[i]);
	}

	for(i=0;i<4;i++)
		solex[i]=(1/(2*M_PI))*log(1.0/sqrt(Xb[i]*Xb[i]+Yb[i]*Yb[i]));

	// We compute sol=-P*neum +M*diri
	
	petscMatVecMult(P,4,N,neumex,Pneumex);
	petscMatVecMult(M,4,N,diriex,Mdiriex);

	for(i=0;i<4;i++)
		sol[i]=-Pneumex[i]+Mdiriex[i];

	double err=0;
	for(i=0;i<4;i++){
		err=err+pow((solex[i]-sol[i])/solex[i],2);
	}
	
	err=err/N;
	err=sqrt(err);
	return err;
}

double test_peanutBulkDx(int N)
{

	// geometric data of the boundary
	double *X,*Y,*xn,*yn,*nx,*ny,*H;


	// points inside the bulk
	
	double Xb[4]={2,2,-2,-2};
	double Yb[4]={2,-2,2,-2};

	int i,j,k,jp;

	createPeanut(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	M=malloc(4*N*sizeof(double));
	P=malloc(4*N*sizeof(double));

	for(i=0;i<4;i++){
		for(j=0;j<N;j++){
			k=j+i*N;
			jp=mod(j+1,N);
			GreenBoundDx(Xb[i],Yb[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],M+k,P+k);
		}
	}

	// Prepare the Dirichlet and Neumann data, compute the exact solution.

	double *neumex,*diriex;
	double *Pneumex,*Mdiriex;
	double *solexDx,*solDx;

	Mdiriex=malloc(4*sizeof(double));
	diriex=malloc(N*sizeof(double));
	Pneumex=malloc(4*sizeof(double));
	neumex=malloc(N*sizeof(double));
	solexDx=malloc(4*sizeof(double));
	solDx=malloc(4*sizeof(double));

	for(i=0;i<N;i++){
		diriex[i]=(1/(2*M_PI))*log(1.0/sqrt(xn[i]*xn[i]+yn[i]*yn[i]));
		neumex[i]=-(1/(2*M_PI))*(1.0/(xn[i]*xn[i]+yn[i]*yn[i]))*(xn[i]*nx[i]+yn[i]*ny[i]);
	}

	for(i=0;i<4;i++)
		solexDx[i]=-(1/(2*M_PI))*Xb[i]/(Xb[i]*Xb[i]+Yb[i]*Yb[i]);

	// We compute sol=-P*neum +M*diri
	
	petscMatVecMult(P,4,N,neumex,Pneumex);
	petscMatVecMult(M,4,N,diriex,Mdiriex);

	for(i=0;i<4;i++)
		solDx[i]=-Pneumex[i]+Mdiriex[i];

	double err=0;
	for(i=0;i<4;i++){
		err=err+pow((solexDx[i]-solDx[i])/solexDx[i],2);
	}
	
	err=err/N;
	err=sqrt(err);
	return err;
}


double test_peanutBulkDy(int N)
{

	// geometric data of the boundary
	double *X,*Y,*xn,*yn,*nx,*ny,*H;


	// points inside the bulk
	
	double Xb[4]={2,2,-2,-2};
	double Yb[4]={2,-2,2,-2};

	int i,j,k,jp;

	createPeanut(N,&X,&Y,&xn,&yn,&nx,&ny,&H);

	// Here, we build the influence matrices M and P

	double* M;
	double* P;

	M=malloc(4*N*sizeof(double));
	P=malloc(4*N*sizeof(double));

	for(i=0;i<4;i++){
		for(j=0;j<N;j++){
			k=j+i*N;
			jp=mod(j+1,N);
			GreenBoundDy(Xb[i],Yb[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],M+k,P+k);
		}
	}

	// Prepare the Dirichlet and Neumann data, compute the exact solution.

	double *neumex,*diriex;
	double *Pneumex,*Mdiriex;
	double *solexDy,*solDy;

	Mdiriex=malloc(4*sizeof(double));
	diriex=malloc(N*sizeof(double));
	Pneumex=malloc(4*sizeof(double));
	neumex=malloc(N*sizeof(double));
	solexDy=malloc(4*sizeof(double));
	solDy=malloc(4*sizeof(double));

	for(i=0;i<N;i++){
		diriex[i]=(1/(2*M_PI))*log(1.0/sqrt(xn[i]*xn[i]+yn[i]*yn[i]));
		neumex[i]=-(1/(2*M_PI))*(1.0/(xn[i]*xn[i]+yn[i]*yn[i]))*(xn[i]*nx[i]+yn[i]*ny[i]);
	}

	for(i=0;i<4;i++)
		solexDy[i]=-(1/(2*M_PI))*Yb[i]/(Xb[i]*Xb[i]+Yb[i]*Yb[i]);

	// We compute sol=-P*neum +M*diri
	
	petscMatVecMult(P,4,N,neumex,Pneumex);
	petscMatVecMult(M,4,N,diriex,Mdiriex);

	for(i=0;i<4;i++)
		solDy[i]=-Pneumex[i]+Mdiriex[i];

	double err=0;
	for(i=0;i<4;i++){
		err=err+pow((solexDy[i]-solDy[i])/solexDy[i],2);
	}
	
	err=err/N;
	err=sqrt(err);
	return err;
}
// Running all the test and print the info

int main(){
	int i;
	petscInit();
	int N[4]={50,100,200,400};

	for(i=0;i<4;i++)
		printf("test_peanutDN:\t n=%i\terr=%e\n",N[i],test_peanutDN(N[i]));
	for(i=0;i<4;i++)
		printf("test_peanutND:\t n=%i\terr=%e\n",N[i],test_peanutND(N[i]));

	// tests in the bulk

	int M[5]={50,100,200,400,800};

	for(i=0;i<5;i++)
		printf("test_peanutBulk:\t n=%i\terr(relative)=%e\n",M[i],test_peanutBulk(M[i]));
	for(i=0;i<5;i++)
		printf("test_peanutBulkDx:\t n=%i\terr(relative)=%e\n",M[i],test_peanutBulkDx(M[i]));
	for(i=0;i<5;i++)
                printf("test_peanutBulkDy:\t n=%i\terr(relative)=%e\n",M[i],test_peanutBulkDy(M[i]));
	petscEnd();
	return 0;
}

