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
