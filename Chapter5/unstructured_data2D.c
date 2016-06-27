
#include <stdlib.h>
#include <stdio.h>

#ifndef N
#define N 32 
#endif
#ifndef M
#define M 32 
#endif

double ** allocData(size_t size1, size_t size2);
int deleteData(double ** A, size_t size1);
int initData(double **A, size_t size1, size_t size2, double val);
int printData(double **A, size_t size1, size_t size2);

int main() {

    double **A, **B;
    size_t size1,size2, i, j;
    size1 = N;
    size2 = M;
    A=allocData(size1,size2);
    B=allocData(size1,size2);
    initData(B,size1,size2,2.5);

/* Perform the computation on the device */
#pragma acc parallel loop collapse(2) present(A,B)
    for (j=0; j < size1; ++j) {
       for (i=0; i < size2; ++i) {
          A[j][i] = B[j][i] + (double) ((j*size2)+i);
       }
    } 
/* Copy back the results */ 
#pragma acc update self(A[0:size1][0:size2])
    printData(A,size1,size2);
    deleteData(A,size1);
    deleteData(B,size1);
    exit(0);
}

double ** allocData(size_t size1, size_t size2) {
    double ** tmp;
    int i;
    tmp = (double **) malloc(size1*sizeof(double*));
    for (i=0; i < size1; ++i) {
       tmp[i] = (double *) malloc(size2*sizeof(double)); 
    }
#pragma acc enter data create(tmp[0:size1][0:size2])
    return tmp;
}

int deleteData(double ** A, size_t size1) {
    int i;
#pragma acc exit data delete(A)
    for (i=0; i < size1; ++i) {
       free(A[i]);
    }
    free(A);
}

int initData(double **A, size_t size1, size_t size2, double val) {
    size_t i,j;
    for (j=0; j < size1; ++j) {
       for (i=0; i < size2; ++i) {
	  A[j][i] = val;
       }
    }
/* Update the device with the initial values */
#pragma acc update device(A[0:size1][0:size2])
}

int printData(double **A, size_t size1, size_t size2) {
    size_t i,j,n;
    n = size2-1;
    printf("Values:\n");
    for (i=0; i < 5; ++i) {   
	printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i][0],i,n,A[i][n]);
    } 
    printf("....\n");
    for (i=size1-5; i < size1; ++i) {   
	printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i][0],i,n,A[i][n]);
    } 
}


