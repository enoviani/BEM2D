#include"petscInterface.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Just to print out the result of the tests

const char* isPass(int i){
        static char str1[]="pass";
	static char str2[]="no pass";
	if(i) 
		return str1;
	else 
		return str2;
}

// A simple test for linear equation solving

int test_petscSolve(){
	double M[4]={2,-1,-1,2};
	double rhs[2]={1,1};
	double x[2]={0,0};
	petscSolve(M,2,rhs,x);
	int pass;
	if(sqrt((x[0]-1)*(x[0]-1)+(x[1]-1)*(x[1]-1))<1e-5)
		pass=1;
	else
		pass=0;
	return pass;
}

int test_petscMatVecMult(){
	double M[6]={1,2,3,4,5,6};
	double X[3]={1,0,1};
	double Y[2]={0,0};
	petscMatVecMult(M,2,3,X,Y);
	int pass;
	if(sqrt((Y[0]-4)*(Y[0]-4)+(Y[1]-10)*(Y[1]-10))<1e-5)
		pass=1;
	else
		pass=0;
	return pass;
}

// Running all the test and print the pass/nopass statement

int main(){
	petscInit();
	printf("test_petscSolve:\t%s\n",isPass(test_petscSolve()));
	printf("test_petscMatVecMult:\t%s\n",isPass(test_petscMatVecMult()));
	petscEnd();
	return 0;
}

