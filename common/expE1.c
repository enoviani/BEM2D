#include "expE1.h"

double complex ExpE1(double complex z){
	int i;
	double complex ExE1;
	double Re, Ima;
	double complex overmZ;
	double complex Temp;
	double nu, m1,m2,m3,m4,n1,n2,n3,n4,n5,d1,d2,d3,d4,d5,d6;
	double complex M;
	double complex N;
	double complex D;

	//	printf("%f,%f",creal(M),cimag(M));
	Re=creal(z);
	Ima=cimag(z);
	if(Re<=-16){
		overmZ=1/(-z);
		Temp=1/z;
		ExE1=Temp;
		for(i=1;i<6;i++){
			Temp=Temp*i*overmZ;
			ExE1=ExE1+Temp;
		}
	}


	else if(Re>-0.5){
		nu=0.5772156649;
		m1=0.23721365;
		m2=0.0206543;
		m3=0.000763297;
		m4=0.0000097087007;

		n1=-1.49545886;
		n2=0.041806426;
		n3=-0.03000591;
		n4=0.0019387339;
		n5=-0.00051801555;

		d1=-0.76273617;
		d2=0.28388363;
		d3=-0.066786033;
		d4=0.012982719;
		d5=-0.0008700861;
		d6=0.0002989204;


		M=-(1+m1*z+m2*cpow(z,2)+m3*cpow(z,3)+m4*cpow(z,4))*clog(z);
		N=-nu*(0.99999207+n1*z+n2*cpow(z,2)+n3*cpow(z,3)+n4*cpow(z,4)+n5*cpow(z,5));
		D=1+d1*z+d2*cpow(z,2)+d3*cpow(z,3)+d4*cpow(z,4)+d5*cpow(z,5)+d6*cpow(z,6);
		//	printf("%f,%f",creal(M),cimag(M));
		ExE1=(M+N)/D;
	}
	else if(Re>10 || abs(Ima)>10){

		ExE1=0.711093/(z+0.415775) + 0.278518/(z+2.29428) + 0.010389/(z+6.29);
	}
	else{
		// Zhang's formula instead of Delhommeau
		//
		double el = 0.5772156649015328; 
		double x = creal(z);
		double a0 = cabs(z);
		double complex ce1;

		if(a0==0.0)
			ce1 = 1e300; 
		else if( a0 <= 10.0 || ( x < 0 && a0 < 20) ){
			ce1 = 1;
			double complex cr = 1;
			int k;
			for(k=1;k<=150;k++){
				cr = -cr*k*z/cpow(k+1.0,2);
				ce1 = ce1 + cr;
				if ( cabs (cr) <= cabs(ce1)*1.0e-15)
					break;
			}
			ce1=-el-clog(z)+z*ce1;

		}
		else{
			double complex ct0 = 0;
			int k;
			for(k=120;k>=1;k--)
				ct0=k/(1+k/(z+ct0));
			double complex ct=1.0/(z+ct0);
			ce1=cexp(-z)*ct;
			if ( x <= 0 && cimag(z)==0)
				ce1 = ce1 - I*M_PI;
		}
		ExE1=ce1*cexp(z);
	}
	return ExE1;
}

