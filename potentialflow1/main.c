#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"green.h"
#include"display.h"
#include<complex.h>
#include"emptyMatrix.h"
#include"matrixOperation.h"
#include"petscInterface.h"
#include <petscksp.h>


int mod (int a, int b)
{
	int ret = a % b;
	if(ret < 0)
		ret+=b;
	return ret;
}



int main(){

	// Make the Obstacle, a circle
	double* X;
	double* Y;
	double* sudut;
	double PI=3.14159265;
	double Npart=10; //number of partition 
	int N=Npart;
	double dt=1.0/((double)N);;
	int i,j,ip,jp;
	double** IdPlusM;

	X=malloc(N*sizeof(double));
	Y=malloc(N*sizeof(double));
	sudut=malloc(N*sizeof(double));
	dt=1/Npart;

	for(i=0;i<N;i++){
		sudut[i] = 2*PI*i*dt;
		X[i]= cos(sudut[i]);
		Y[i]= sin(sudut[i]);
	}

	printf("Coordinates of Obstacle \n");
	displayCoordinates(X,Y,N);
	//displayArray(Y,N);

	// xi
	double* xn;
	double* yn;
	xn = malloc(N*sizeof(double));
	yn = malloc(N*sizeof(double));

	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		xn[i]= (X[ip]+X[i])/2;
		yn[i]= (Y[ip]+Y[i])/2;
	}

	printf("coordinates of Xi\n");
	displayCoordinates(xn,yn,N);
	// displayArray(yn,N);

	//normal vector
	double* nx;
	double* ny;
	double* norm;
	nx=malloc(N*sizeof(double));
	ny=malloc(N*sizeof(double));
	norm=malloc(N*sizeof(double));

	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		nx[i]= -(Y[ip]-Y[i]);
		ny[i]= X[ip]-X[i];
		norm[i]= sqrt(pow(nx[i],2)+pow(ny[i],2));
		nx[i]=nx[i]/norm[i];
		ny[i]=ny[i]/norm[i];
	}

	printf("normal vector \n");
	displayCoordinates(nx,ny,N);
	printf("xi vector and normal vector \n");
	displayNormal(xn,yn,nx,ny,N);


	// Compute matrices M and P
	double** M;
	double** P;
	double* li;

	li= malloc(N*sizeof(double));
	M=createEmptyMatrix(N,N);
	P=createEmptyMatrix(N,N);


	for(i=0;i<N;i++){
		ip=mod(i+1,N);
		li[i]= sqrt(pow((X[ip]-X[i]),2)+pow((Y[ip]-Y[i]),2));
	}
	printf("li: \n");
	displayArray(li,N);

	//Create Identity matrix I (the matrix needed is 0.5*I)
	double** Id;
	printf("matrix I is \n");
	Id= createEmptyMatrix(N,N);

	for(i=0;i<N;i++){
		for (j=0;j<N;j++){
			if (i==j){
				Id[i][j]=0.5;
			}
			else { 
				Id[i][j]=0.0;
			}
		}
	}

	displayMatrix(Id,N);  

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			jp=mod(j+1,N);
			if (i==j){
				M[i][j]=0.0;
				P[i][j]=(1/(2*PI))*((li[i]*log(li[i]/2))-li[i]);
			}
			else{
				GreenBound(xn[i],yn[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],&(M[i][j]),&(P[i][j]));
			}
		}

	}
	printf("Matrix M is:\n");
	displayMatrix(M,N);
	printf("Matrix P is:\n");
	displayMatrix(P,N);


	//The next step is finding the solution of linear equation.

	//Psi=(0.5*I+M) \ P*(nx');


	IdPlusM=createEmptyMatrix(N,N);
	addMatrix(Id, N, N, M,N, N, IdPlusM);

	double* Pn;
	Pn=malloc(N*sizeof(double));
	multipleMatrix(P, N, N, nx, N, Pn);	
	printf("we will solve (0.5*I+M)*Psi= P*(nx')\n");
	printf("RHS :\n");
	displayArray(Pn,N);
	printf("The result of (0.5*I+M): \n");
	displayMatrix(IdPlusM,N); 

	//We will solveIdPlusM*Psi=Pn

	double* A;
	A= malloc((N*N)*sizeof(double));
	//Change matrix (double**) to array (double*)
	int k;
	k=0;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[k]=IdPlusM[i][j];
			k++;
		}
	}
	printf("matrix IdPlusM in lines order\n");
	displayArray(A,(N*N));
	double* x1;
	x1= malloc(N*sizeof(double));
	petscSolve(A, N, Pn, x1);//Solving the Linear Equation System
	printf("Solution Psi:\n");
	displayArray(x1,N);

	//**************************************************
	//      Calculate potential in all the Domain
	//**************************************************
	/*
	//Building Discretization of domain

	int Hx=5;
	int Hy=3;

	double dhx=(1-(-1))/(double)Hx; //Domain x is in [-1,1] with distance between successor and predesessor at x is dhx
	double dhy=(1-(-1))/(double)Hy; //Domain y is in [-1,1] with distance between successor and predesessor at y is dhy

	double* XXV;
	double* YYV;

	XXV=malloc((Hx*Hy)*sizeof(double));
	YYV=malloc((Hx*Hy)*sizeof(double));

	int k=0;

	for (i=0;i<Hy;i++){
	k=i*Hx;
	for(j=0;j<Hx;j++){
	XXV[k+j]=(-1)+dhx*j;
	YYV[k+j]=(-1)+dhy*i;
	}
	}
	printf("vector XXV :\n");
	displayArray(XXV,(Hx*Hy));

	printf("vector YYV :\n");
	displayArray(YYV,(Hx*Hy));



	 */

	return 1;

}
