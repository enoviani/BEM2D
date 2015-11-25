#include"green.h"


void GreenBoundDiag(double H, double* pm, double* pp){
	*pm=0.5;
	*pp=(1/(2*M_PI))*(H*log(H/2)-H);
}

void GreenBound(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pm, double* pp){

  double amx=xm-xx;
  double amy=ym-yy;

  double normam=sqrt(pow(amx,2)+pow(amy,2));

  amx=amx/normam;
  amy=amy/normam;

  double apx=xp-xx;
  double apy=yp-yy;

  double normap=sqrt(pow(apx,2)+pow(apy,2));

  apx=apx/normap;
  apy=apy/normap;

  double sig=copysign(1.0,(apx*nx+apy*ny));

  double ThetaM = copysign(1.0,(sig*(nx*amy-ny*amx)))*acos(sig*(nx*amx+ny*amy));
  double ThetaP = copysign(1.0,(sig*(nx*apy-ny*apx)))*acos(sig*(nx*apx+ny*apy));

  double  H = cos(ThetaM)*sqrt(pow(xm-xx,2)+pow(ym-yy,2));

//  printf("%f \t %f \n",copysign(1.0,(sig*(nx*amy-ny*amx))),copysign(1.0,(sig*(nx*apy-ny*apx))));

  *pm=(1/(2*M_PI))*(ThetaP-ThetaM);
  *pp=(-sig*(H/(2*M_PI))*((tan(ThetaP)*(log(H*(1.0/cos(ThetaP)))-1)+ThetaP)-(tan(ThetaM)*(log(H*(1.0/cos(ThetaM)))-1)+ThetaM)));
  
}

void GreenBoundDx(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pm, double* pp){

  double amx=xm-xx;
  double amy=ym-yy;

  double normam=sqrt(pow(amx,2)+pow(amy,2));

  amx=amx/normam;
  amy=amy/normam;

  double apx=xp-xx;
  double apy=yp-yy;

  double normap=sqrt(pow(apx,2)+pow(apy,2));

  apx=apx/normap;
  apy=apy/normap;

  double sig=copysign(1.0,(apx*nx+apy*ny));

  double ThetaM = copysign(1.0,(sig*(nx*amy-ny*amx)))*acos(sig*(nx*amx+ny*amy));
  double ThetaP = copysign(1.0,(sig*(nx*apy-ny*apx)))*acos(sig*(nx*apx+ny*apy));

  double  H = cos(ThetaM)*sqrt(pow(xm-xx,2)+pow(ym-yy,2));

//  printf("%f \t %f \n",copysign(1.0,(sig*(nx*amy-ny*amx))),copysign(1.0,(sig*(nx*apy-ny*apx))));


  *pp=(1/(2*M_PI))*(nx*(ThetaP-ThetaM)+ny*(log(fabs(cos(ThetaP)))-log(fabs(cos(ThetaM)))));

  *pm=(1/(2*M_PI*H))*sig*(-nx*0.5*((sin(2*ThetaP)-sin(2*ThetaM)))+ny*(sin(ThetaP)*sin(ThetaP)-sin(ThetaM)*sin(ThetaM))) ; /*new*/
}

void GreenBoundDy(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pm, double* pp){

  double amx=xm-xx;
  double amy=ym-yy;

  double normam=sqrt(pow(amx,2)+pow(amy,2));

  amx=amx/normam;
  amy=amy/normam;

  double apx=xp-xx;
  double apy=yp-yy;

  double normap=sqrt(pow(apx,2)+pow(apy,2));

  apx=apx/normap;
  apy=apy/normap;

  double sig=copysign(1.0,(apx*nx+apy*ny));

  double ThetaM = copysign(1.0,(sig*(nx*amy-ny*amx)))*acos(sig*(nx*amx+ny*amy));
  double ThetaP = copysign(1.0,(sig*(nx*apy-ny*apx)))*acos(sig*(nx*apx+ny*apy));

  double  H = cos(ThetaM)*sqrt(pow(xm-xx,2)+pow(ym-yy,2));

  *pp=(1/(2*M_PI))*(ny*(ThetaP-ThetaM)-nx*(log(fabs(cos(ThetaP)))-log(fabs(cos(ThetaM)))));
  *pm=sig*(1/(2*M_PI*H))*(-ny*0.5*((sin(2*ThetaP)-sin(2*ThetaM)))-nx*(sin(ThetaP)*sin(ThetaP)-sin(ThetaM)*sin(ThetaM)));

 // printf("%e \t %e\n",cos(ThetaP),cos(ThetaM));
}

void GreenRegularH(double xx, double yy, double xm, double ym, double xp,double yp, double nx, double ny,double k, double* mDH, double* pH){
	double rew;
	double imw;
	double L;

	rew = -(0.5*(ym+yp)+yy);
	imw = 0.5*(xm+xp)-xx;
	//printf("w=%f\t+\t%fi\n",rew,imw);

	double complex w = rew + imw * I;
	
	L=sqrt(pow((xp-xm),2)+pow((yp-ym),2));

	if(imw<0){
		*pH=(-(2/M_PI)*creal(ExpE1(-k*w) + clog(-k*w)))*L; 
		*mDH=((-2*k/M_PI)*cimag(ExpE1(-k*w))*nx+(-2*k/M_PI)*creal(ExpE1(-k*w))*ny)*L; 
	}
	else{
		*pH=(-(2/M_PI)*creal(ExpE1(-k*w) +cexp(-k*w)*2*M_PI*I + clog(-k*w)))*L; 
		*mDH=((-2*k/M_PI)*cimag(ExpE1(-k*w)+cexp(-k*w)*2*M_PI*I)*nx+(-2*k/M_PI)*creal(ExpE1(-k*w)+cexp(-k*w)*2*M_PI*I)*ny)*L; 
	}
}


/*
void GreenRegularHDx(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pmDH, double* ppH){

	double rew;
	double imw;
	double k=1;
	double CE1[2], W[2]; 
	double L;

	rew = -(ym-yy);
	imw = xm-xx;
	//printf("w=%f\t+\t%fi\n",rew,imw);

	double complex w = rew + imw * I;
	double complex overw;

	W[0]=-k*rew; W[1]=-k*imw;
	overw=(1/(pow(rew,2)+pow(imw,2)))*rew+(1/(pow(rew,2)+pow(imw,2)))*imw*I;
	E1Z(W,CE1);
	double complex E1kw = CE1[0] + CE1[1] * I;
	L=sqrt(pow((xp-xm),2)+pow((yp-ym),2));

	if(imw<0){
		*ppH=(2*k/M_PI*cimag((cexp(-k*w)*E1kw)))*L; 
		*pmDH=(-2*k/M_PI)*((nx*(creal((k*cexp(-k*w)*E1kw)+overw)))-(ny*(cimag((k*cexp(-k*w)*E1kw)+overw))))*L;
	}
	else{

		*ppH=(2*k/M_PI*cimag((cexp(-k*w)*(E1kw+2*M_PI*I))))*L;
		*pmDH=(-2*k/M_PI)*((nx*(creal((k*cexp(-k*w)*(E1kw+(2*M_PI*I)))+overw)))-(ny*(cimag((k*cexp(-k*w)*(E1kw+(2*M_PI*I)))+overw))))*L;
	}
}

void GreenRegularHDy(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pmDH, double* ppH){

	double rew;
	double imw;
	double k=1;
	double CE1[2], W[2]; 
	double L;

	rew = -(ym-yy);
	imw = xm-xx;
	//printf("w=%f\t+\t%fi\n",rew,imw);

	double complex w = rew + imw * I;
	double complex overw;

	W[0]=-k*rew; W[1]=-k*imw;
	overw=(1/(pow(rew,2)+pow(imw,2)))*rew+(1/(pow(rew,2)+pow(imw,2)))*imw*I;
	E1Z(W,CE1);
	double complex E1kw = CE1[0] + CE1[1] * I;
	L=sqrt(pow((xp-xm),2)+pow((yp-ym),2));


	if(imw<0){
		*ppH=(-2*k/M_PI*creal((cexp(-k*w)*E1kw)))*L; 
		*pmDH=(-2*k/M_PI)*((nx*(cimag((k*cexp(-k*w)*E1kw)+overw)))+(ny*(creal((k*cexp(-k*w)*E1kw)+overw))))*L;
	}
	else{
		*ppH=(-2*k/M_PI*cimag((cexp(-k*w)*(E1kw+2*M_PI*I))))*L;
		*pmDH=(-2*k/M_PI)*((nx*(cimag((k*cexp(-k*w)*(E1kw+(2*M_PI*I)))+overw)))+(ny*(creal((k*cexp(-k*w)*(E1kw+(2*M_PI*I)))+overw))))*L;
	} 

}*/
