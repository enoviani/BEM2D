#include"petscInterface.h"
#include<stdio.h>
#include<stdlib.h>
int main(){
	printf("Hi there ! \n");
	double b[2]={1,1};
	double A[4]={2,1,-1,2};
	double x[2]={0,0};

	petscSolve(A,2,b,x);

	printf("The result is : (%f,%f)\n",x[0],x[1]);

	return 0;
}
