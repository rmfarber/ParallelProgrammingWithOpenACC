#include <stdlib.h>
#include <stdio.h>

#ifndef N
#define N 32 
#endif

double ** allocData(size_t size1, int * sizes);
int deleteData(double ** A, size_t size1);
int initData(double **A, size_t size1, int * sizes, double val);
int printData(double **A, size_t size1, int * sizes);

int main() {

    double **A, **B;
    int * sizes;
    size_t size1,size2, i, j;
    size1 = N;
   
    sizes=(int*) malloc(sizeof(int)*size1); 
    for (i=0; i < size1; ++i) {
       sizes[i]=i+10;
    }
#pragma acc enter data copyin(sizes[0:size1])
    A=allocData(size1,sizes);
    B=allocData(size1,sizes);
    initData(B,size1,sizes,2.5);

/* Perform the computation on the device */
#pragma acc parallel loop gang present(A,B,sizes)
    for (j=0; j < size1; ++j) {
       int size2 = sizes[j]; 
#pragma acc loop vector
       for (i=0; i < size2; ++i) {
          A[j][i] = B[j][i] + (double) ((j*size2)+i);
       }
    } 
#ifdef _OPENACC
/* Copy back the results */ 
    for (j=0; j < size1; ++j) {
      int nele = sizes[j];
#pragma acc update self (A[j:1][0:sizes[j]])
    }
#endif

    printData(A,size1,sizes);
    deleteData(A,size1);
    deleteData(B,size1);
#pragma acc exit data delete(sizes)
    free(sizes);
    exit(0);
}

double ** allocData(size_t size1, int * sizes) {
    double ** tmp;
    int i;
    tmp = (double **) malloc(size1*sizeof(double*));
/* Create an array of pointers */
#pragma acc enter data create(tmp[0:size1][0:1])
    for (i=0; i < size1; ++i) {
       tmp[i] = (double *) malloc(sizes[i]*sizeof(double));
/* Create the vector and attach it to the pointer array */
#pragma acc enter data create(tmp[i:1][0:sizes[i]])
    }
    return tmp;
}

int deleteData(double ** A, size_t size1) {
    int i;
    for (i=0; i < size1; ++i) {
       free(A[i]);
#pragma acc exit data delete(A[i:1])
  }
#pragma acc exit data delete(A)
    free(A);
}

int initData(double **A, size_t size1, int * sizes, double val) {
    size_t i,j;
    for (j=0; j < size1; ++j) {
       int size2=sizes[j];
       for (i=0; i < size2; ++i) {
	  A[j][i] = val;
       }
/* Update the device with the initial values */
#pragma acc update device(A[j:1][0:size2])
    }
}

int printData(double **A, size_t size1, int * sizes) {
    size_t i,j;
    printf("Values:\n");
    for (i=0; i < 5; ++i) {   
        int last = sizes[i]-1;
	printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i][0],i,last,A[i][last]);
    } 
    printf("....\n");
    for (i=size1-5; i < size1; ++i) {   
        int last = sizes[i]-1;
	printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i][0],i,last,A[i][last]);
    } 
}


