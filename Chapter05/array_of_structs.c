#include <stdlib.h>
#include <stdio.h>

#ifndef N
#define N 32 
#endif

typedef struct {
   int size; 
   double * data;
} vector;

vector * allocData(size_t size1);
int deleteData(vector* A, size_t size1);
int initData(vector *A, size_t size1, double val);
int printData(vector *A, size_t size1);

int main() {

    vector *A, *B;
    size_t size1, i, j;
    size1 = N;
   
    A=allocData(size1);
    B=allocData(size1);
    initData(B,size1,2.5);

/* Perform the computation on the device */
#pragma acc parallel loop gang present(A,B)
    for (j=0; j < size1; ++j) {
       int size2 = A[j].size;
#pragma acc loop vector
       for (i=0; i < size2; ++i) {
          A[j].data[i]= B[j].data[i] + (double) ((j*size2)+i);
       }
    } 
#ifdef _OPENACC
/* Copy back the results */ 
    for (j=0; j < size1; ++j) {
#pragma acc update self (A[j].data[0:A[j].size])
    }
#endif

    printData(A,size1);
    deleteData(A,size1);
    deleteData(B,size1);
    exit(0);
}

vector * allocData(size_t size1) {
    vector * tmp;
    int i;
    tmp = (vector*) malloc(size1*sizeof(vector));
/* Create an array of pointers */
#pragma acc enter data create(tmp[0:size1])
    for (i=0; i < size1; ++i) {
       tmp[i].size = i+10;
       tmp[i].data = (double *) malloc(tmp[i].size*sizeof(double));
/* Create the vector and attach it to the pointer array */
#pragma acc enter data create(tmp[i].data[0:tmp[i].size])
/* Update the device's size */
#pragma acc update device(tmp[i].size)
    }
    return tmp;
}

int deleteData(vector * A, size_t size1) {
    int i;
    for (i=0; i < size1; ++i) {
       free(A[i].data);
#pragma acc exit data delete(A[i].data)
  }
#pragma acc exit data delete(A)
    free(A);
}

int initData(vector *A, size_t size1, double val) {
    size_t i,j;
    for (j=0; j < size1; ++j) {
       int size2=A[j].size;
       for (i=0; i < size2; ++i) {
	  A[j].data[i] = val;
       }
/* Update the device with the initial values */
#pragma acc update device(A[j].data[0:size2])
    }
}

int printData(vector *A, size_t size1) {
    size_t i,j;
    printf("Values:\n");
    for (i=0; i < 5; ++i) {   
        int last = A[i].size-1;
	printf("A[%d].data[0]=%f A[%d].data[%d]=%f\n",i,A[i].data[0],i,last,A[i].data[last]);
    } 
    printf("....\n");
    for (i=size1-5; i < size1; ++i) {   
        int last = A[i].size-1;
	printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i].data[0],i,last,A[i].data[last]);
    } 
}


