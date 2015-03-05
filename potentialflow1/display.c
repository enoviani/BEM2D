#include"display.h"

void displayArray(double* array1, int size){
  int i;  
  for(i=0;i<size;i++){
    printf("%f\n",array1[i]);
  }
  printf("----------------\n");
}


void displayCoordinates(double* array1, double* array2, int size){
  int i;  
  for(i=0;i<size;i++){
    
    printf("%f\t%f\n",array1[i],array2[i]);
  }
printf("----------------\n");
}
void displayNormal(double* array1, double* array2,double* array3, double* array4, int size){
  int i;  
  for(i=0;i<size;i++){
    
    printf("%f\t%f\t%f\t%f\n",array1[i],array2[i],array3[i],array4[i]);
  }
printf("----------------\n");
}


void displayMatrix(double** matrix, int N){
  int i,j;
    printf("----------\n");
  for(i=0;i<N;i++){
    
    for(j=0;j<N;j++){
      printf("%f\t",matrix[i][j]);
    }
    printf("\n");
  }
}
