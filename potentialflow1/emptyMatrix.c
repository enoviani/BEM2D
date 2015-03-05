#include"emptyMatrix.h"


double** createEmptyMatrix(int M,int N){
  double** Matrix= malloc(M*sizeof(double*));
  int i;
  for(i=0;i<M;i++){
    Matrix[i]=malloc(N*sizeof(double));
	
  }
  return Matrix;
}
