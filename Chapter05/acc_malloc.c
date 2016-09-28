#include <stdlib.h>
#include <stdio.h>

#ifdef _OPENACC
#include <openacc.h>
#endif

#ifndef N
#define N 1024
#endif

int printData(double *, size_t ); 
void cudaSet(double *, size_t, double);

int main() {
    double *A, *B;
    size_t size, i;
    size = N;
    A= (double*) malloc(size*sizeof(double));
    B= (double*) acc_malloc(size*sizeof(double));

    /* Call a CUDA routine to set the device data */
    cudaSet(B,size,2.5);

    /* Use the deviceptr clause on the compute region */
    #pragma acc parallel loop copyout(A[0:size]) deviceptr(B)
    for (i=0; i < size; ++i) {
       A[i] = B[i] + (double) i;
    }
    printData(A, size);
    free(A);  
    acc_free(B); 
    exit(0);
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

