#include"green.h"


void GreenBoundDiag(double H, double* pm, double* pp){
	*pm=-0.5;
	*pp=(1/(2*M_PI))*(H*log(H/2)-H);
}

void GreenBound(double xx, double yy, double xm, double ym,double xp,double yp, double nx, double ny, double* pm, double* pp){

  double PI=3.14159265;
  
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

  *pm=(1/(2*PI))*(ThetaP-ThetaM);
  *pp=(-sig*(H/(2*PI))*((tan(ThetaP)*(clog(H*(1.0/cos(ThetaP)))-1)+ThetaP)-(tan(ThetaM)*(clog(H*(1.0/cos(ThetaM)))-1)+ThetaM)));
  
}
