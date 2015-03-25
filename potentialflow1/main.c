#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"green.h"
#include"display.h"
#include<complex.h>
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
	double Npart=100; //number of partition 
	int N=Npart;
	double dt=1.0/((double)(N-1));;
	int i,j,ip,jp;
	double** IdPlusM;
	double a=1.0; // radius of circle(obstacle)
	double u=2;// velocity far from the obstacle

	petscInit();

	X=malloc(N*sizeof(double));
	Y=malloc(N*sizeof(double));
	sudut=malloc(N*sizeof(double));
	dt=1/Npart;

	for(i=0;i<N;i++){
		sudut[i] = 2*PI*i*dt;
		X[i]=a*(cos(sudut[i]));
		Y[i]=a*(sin(sudut[i]));
	}

//	printf("Coordinates of Obstacle \n");
//	displayCoordinates(X,Y,N);
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

//	printf("coordinates of Xi\n");
//	displayCoordinates(xn,yn,N);
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

//	printf("normal vector \n");
//	displayCoordinates(nx,ny,N);
//	printf("xi vector and normal vector \n");
//	displayNormal(xn,yn,nx,ny,N);


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
/*	printf("li: \n");
	displayArray(li,N);
*/
	//Create Identity matrix I (the matrix needed is 0.5*I)
	double** Id;
//	printf("matrix I is \n");
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

//	displayMatrix(Id,N);  

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
/*	printf("Matrix M is:\n");
	displayMatrix(M,N);
	printf("Matrix P is:\n");
	displayMatrix(P,N);

*/
	//The next step is finding the solution of linear equation.

	//Psi=(0.5*I+M) \ P*(nx');


	IdPlusM=createEmptyMatrix(N,N);
	addMatrix(Id, N, N, M,N, N, IdPlusM);

	double* Pn;
	Pn=malloc(N*sizeof(double));
	multipleMatrix(P, N, N, nx, N, Pn);	
//	printf("we will solve (0.5*I+M)*Psi= P*(nx')\n");
//	printf("RHS :\n");
//	displayArray(Pn,N);
//	printf("The result of (0.5*I+M): \n");
//	displayMatrix(IdPlusM,N); 

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
//	printf("matrix IdPlusM in lines order\n");
//	displayArray(A,(N*N));
	double* Psi;
	Psi= malloc(N*sizeof(double));
	petscSolve(A, N, Pn, Psi);//Solving the Linear Equation System
//	printf("Solution Psi:\n");
//	displayArray(Psi,N);

	//**************************************************
	//      Calculate potential in all the Domain
	//**************************************************

	//Building Discretization of domain

	int Hx=11;
	int Hy=11;
	int numpoint=Hx*Hy;

	double mesh1=-1.2;
	double mesh2=1.2;
	double dhx=(mesh2-mesh1)/(double)(Hx-1); //Domain x is in [-1,1] with distance between successor and predesessor at x is dhx
	double dhy=(mesh2-mesh1)/(double)(Hy-1); //Domain y is in [-1,1] with distance between successor and predesessor at y is dhy

	double* XXV;
	double* YYV;
	// Building mesh in array (double*) //
	XXV=malloc((numpoint)*sizeof(double));
	YYV=malloc((numpoint)*sizeof(double));

	k=0;
	for (i=0;i<Hx;i++){
		k=i*Hx;
		for(j=0;j<Hy;j++){
			XXV[k+j]=(mesh1)+dhx*i;
			YYV[k+j]=(mesh1)+dhy*j;
		}
	}
//	printf("vector XXV :\n");
//	displayArray(XXV,(Hx*Hy));

//	printf("vector YYV :\n");
//	displayArray(YYV,(Hx*Hy));

	double** MMx;
	double** PPx;
	double** MMy;
	double** PPy;
	double** MM;
	double** PP;
	MMx=createEmptyMatrix(numpoint,N);
	PPx=createEmptyMatrix(numpoint,N);

	MMy=createEmptyMatrix(numpoint,N);
	PPy=createEmptyMatrix(numpoint,N);

	MM=createEmptyMatrix(numpoint,N);
	PP=createEmptyMatrix(numpoint,N);

	double eps=1E-6L;
	for(i=0;i<numpoint;i++){
		for (j=0;j<N;j++){
			jp=mod(j+1,N);
			GreenBound(XXV[i],YYV[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],&(MM[i][j]),&(PP[i][j]));

			GreenBound(XXV[i]+eps,YYV[i],X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],&(MMx[i][j]),&(PPx[i][j]));

			GreenBound(XXV[i],YYV[i]+eps,X[j],Y[j],X[jp],Y[jp],nx[j],ny[j],&(MMy[i][j]),&(PPy[i][j]));
			MMx[i][j]=(MMx[i][j]-MM[i][j])/eps;
			PPx[i][j]=(PPx[i][j]-PP[i][j])/eps;
			MMy[i][j]=(MMy[i][j]-MM[i][j])/eps;
			PPy[i][j]=(PPy[i][j]-PP[i][j])/eps;
		}
	}
	//Calculate Psi in domain: int Green*normalx -Psi*int Del Green
	double* PPn;
	double* MMPsi;
	double* PPxn;
	double* MMxPsi;
	double* PPyn;
	double* MMyPsi;
	PPn=malloc(numpoint*sizeof(double));
	MMPsi=malloc(numpoint*sizeof(double));
	PPxn=malloc(numpoint*sizeof(double));
	MMxPsi=malloc(numpoint*sizeof(double));
	PPyn=malloc(numpoint*sizeof(double));
	MMyPsi=malloc(numpoint*sizeof(double));

	multipleMatrix(PPx, numpoint, N, nx, N, PPxn);
	multipleMatrix(MMx, numpoint, N, Psi, N, MMxPsi);
	multipleMatrix(PPy, numpoint, N, nx, N, PPyn);
	multipleMatrix(MMy, numpoint, N, Psi, N, MMyPsi);
	multipleMatrix(PP, numpoint, N, nx, N, PPn);
	multipleMatrix(MM, numpoint, N, Psi, N, MMPsi);

	double* VVx;
	double* VVy;
	double* PPhiv;
	VVx=malloc(numpoint*sizeof(double));
	VVy=malloc(numpoint*sizeof(double));
	PPhiv=malloc(numpoint*sizeof(double));


	for(i=0;i<numpoint;i++){
		VVx[i]=PPxn[i]-MMxPsi[i];
		VVy[i]=PPyn[i]-MMyPsi[i];
		PPhiv[i]=PPn[i]-MMPsi[i];
	}
	printf("*******************************\n X\tY\tPsi: in boundary.txt\n*******************************\n\n\n");
//	displayBoundary3(X,Y,Psi,N);
//write the result Psi at boundary in boundary.txt
	FILE *fp1;

        fp1 = fopen("boundary.txt", "w");

        for(i=0;i<N;i++){

                fprintf(fp1,"%f\t%f\t%f\n",X[i],Y[i],Psi[i]);
        }

        fclose(fp1);


	printf("*****************************\nvalue in Domain:");
	printf("written in domain.txt\n*****************************\n");
//	displayDomain5(XXV,YYV,VVx,VVy, PPhiv, (Hx*Hy));
	//**** write the result to file domain.txt****
	FILE *fp;

        fp = fopen("domain.txt", "w");

	for(i=0;i<numpoint;i++){

                fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",XXV[i],YYV[i],VVx[i]+u,VVy[i],PPhiv[i]);
        }

	fclose(fp);

	/* ********************************************
	   Validate the result
	 ***********************************************/
	for(i=0;i<numpoint;i++){
		VVx[i]=VVx[i]+u;//the actual velocity potential
	}

	double* valVVx;
	double* valVVy;
	double T;
	double r; // r^2=x^2+y^2
	double* errx;
	double* erry;
	valVVx=malloc(numpoint*sizeof(double));
	valVVy=malloc(numpoint*sizeof(double));
	errx=malloc(numpoint*sizeof(double));
	erry=malloc(numpoint*sizeof(double));

	for(i=0;i<numpoint;i++){
		r=sqrt(pow(XXV[i],2)+pow(YYV[i],2));
		T=atan(YYV[i]/XXV[i]);
/*		valVVx[i]=-minu*(r+((pow(a,2)/r)))*cos(T);
		valVVy[i]=-minu*(r+((pow(a,2)/r)))*cos(T);//masih salah
	
		valVVx[i]=u*((pow(a,2)*(pow(YYV[i],2)-pow(XXV[i],2)))+(pow((pow(XXV[i],2)+pow(YYV[i],2)),2)))/(pow((pow(XXV[i],2)+pow(YYV[i],2)),2));
		valVVy[i]=-u*(2*pow(a,2)*XXV[i]*YYV[i])/(pow((pow(XXV[i],2)+pow(YYV[i],2)),2));

*/	
		valVVx[i]= ((-1/2)*u*pow(a,3)*cos(T))/pow(r,2);
		valVVy[i]=(-1/2)*u*pow(r,2)*pow(sin(T),2)*(1-(pow(a,3)/pow(r,3)));
		errx[i]=(abs((valVVx[i]-VVx[i])/valVVx[i]))*100;
		erry[i]=(abs((valVVy[i]-VVy[i])/valVVy[i]))*100;
		printf("%d\t%f\t%f\n",i, errx[i],erry[i]);
	}
	petscEnd();
	return 1;
}
