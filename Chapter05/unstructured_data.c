
#include <stdlib.h>
#include <stdio.h>

#ifndef N
#define N 1024
#endif

double * allocData(size_t size);
int deleteData(double * A);
int initData(double *A, size_t size, double val);
int printData(double *A, size_t size);

int main() {

    double *A, *B;
    size_t size, i;
    size = N;
    A=allocData(size);
    B=allocData(size);
    initData(B,size,2.5);

/* Perform the computation on the device */
#pragma acc parallel loop present(A,B)
    for (i=0; i < size; ++i) {
       A[i] = B[i] + (double) i;
    }
/* Copy back the results */ 
#pragma acc update self(A[0:size])
    printData(A, size);
    deleteData(A);
    deleteData(B);
    exit(0);
}

double * allocData(size_t size) {
    double * tmp;
    tmp = (double *) malloc(size*sizeof(double));
/* Create the array on device. 
   Order matters.  The host copy must be allocated before 
   creating the device copy   */
#pragma acc enter data create(tmp[0:size])
    return tmp;
}

int deleteData(double * A) {
/* Delete the host copy.
   Order matters.  The device copy must be deleted before
   the host copy is freed.  */
#pragma acc exit data delete(A)
    free(A);
}

int initData(double *A, size_t size, double val) {
    size_t i;
    for (i=0; i < size; ++i) {
	A[i] = val;
    }
/* Update the device with the initial values */
#pragma acc update device(A[0:size])
}

int printData(double *A, size_t size) {
    size_t i;
    printf("Values:\n");
    for (i=0; i < 10; ++i) {   
	printf("A[%d]=%f\n",i,A[i]);
    } 
    printf("....\n");
    for (i=size-10; i < size; ++i) {   
	printf("A[%d]=%f\n",i,A[i]);
    } 
}


