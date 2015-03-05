#include"matrixOperation.h"

void addMatrix(double** matrix1, int row1, int col1, double** matrix2,int row2, int col2, double** matrixAdd){
        int i,j;
        if(row1==row2 && col1==col2){
                for (i=0;i<row1;i++){
                        for(j=0;j<col1;j++){
                                matrixAdd[i][j]=matrix1[i][j]+matrix2[i][j];
                        }
                }
        }
        else{
                printf("the size of matrix does not same!!");
        }
}

void multipleMatrix(double** matrix1, int row1, int col1, double* vector,int row2, double* matrixMultiple){
	int i,j,k;
	if (col1==row2){
		for (i=0;i<row1;i++){
			double sum=0.0;
			for(k=0;k<col1;k++){
				sum = sum + (matrix1[i][k]*vector[k]);
			}
			matrixMultiple[i]=sum;
		}
        }
        else
        {printf("the size of matrix is not compatible for mutiplying matrix!!");}
}
