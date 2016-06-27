#include <stdlib.h>
#include <stdio.h>

#include <openacc.h>

#ifndef N
#define N 1024
#endif

#ifndef BLOCKS
#define BLOCKS 32 
#endif

double * allocData(int size);
void initData(double *, int, int, double);
int printData(double **, int, int);

int main() {

    double **A, *tmpA;
    char * devMem;
    int i, j, nbytes;
    nbytes = N*sizeof(double); 

    A= (double**) malloc(BLOCKS*sizeof(double*));
    for (i=0; i < BLOCKS; ++i) {
       A[i] = (double*) malloc(nbytes);
    }
    devMem= (char*) acc_malloc(nbytes);

    for (i=0; i < BLOCKS; ++i) {
	tmpA = A[i];
	acc_map_data(tmpA,devMem,nbytes);
        #pragma acc parallel loop present(tmpA)
        for (j=0; j < N; ++j) {
             tmpA[j] = 2.5 + (double) ((j*N)+i);
        }
	#pragma acc update host(tmpA[0:N])
 	acc_unmap_data(A[i]);
    }  
	
    printData(A,BLOCKS,N);
    for (i=0; i < BLOCKS; ++i) {
	free(A[i]);
    }
    free(A);
    acc_free(devMem);
    exit(0);
 
}

int printData(double **A, int size1, int size2) {
    int i,j,n;
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

