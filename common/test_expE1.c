#include "common.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
	double complex z;
	double complex s;

	printf("Zone 1 : \n");
	z=-20+1*I;
	s=ExpE1(z);
	printf("z\t= %f + i %f\n",creal(z),cimag(z));
	printf("E1(z)\t= %f + i %f\n",creal(s),cimag(s));
	
	printf("Zone 2 : \n");
	z=5+5*I;
	s=ExpE1(z);
	printf("z\t= %f + i %f\n",creal(z),cimag(z));
	printf("E1(z)\t= %f + i %f\n",creal(s),cimag(s));
	
	printf("Zone 3 : \n");
	z=20+20*I;
	s=ExpE1(z);
	printf("z\t= %f + i %f\n",creal(z),cimag(z));
	printf("E1(z)\t= %f + i %f\n",creal(s),cimag(s));

	printf("Zone 4 : \n");
	z=-1+0.1*I;
	s=ExpE1(z);
	printf("z\t= %f + i %f\n",creal(z),cimag(z));
	printf("E1(z)\t= %f + i %f\n",creal(s),cimag(s));

	printf("Large negative real part : \n");
	z=-100000+1*I;
	s=ExpE1(z);
	printf("z\t= %f + i %f\n",creal(z),cimag(z));
	printf("E1(z)\t= %f + i %f\n",creal(s),cimag(s));
	return 0;
}
