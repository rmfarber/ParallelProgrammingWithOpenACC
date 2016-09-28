
#include <stdlib.h>
#include <stdio.h>

#ifndef N
#define N 1024
#endif

typedef struct {
  float x,y,z;
} float3;

float3 * allocData(size_t size);
int deleteData(float3 * A);
int initData(float3 *A, size_t size, float val);
int printData(float3 *A, size_t size);

int main() {

    float3 *A, *B;
    size_t size, i;
    size = N;
    A=allocData(size);
    B=allocData(size);
    initData(B,size,2.5);

/* Perform the computation on the device */
#pragma acc parallel loop present(A,B)
    for (i=0; i < size; ++i) {
       A[i].x = B[i].x + (float) i;
       A[i].y = B[i].y + (float) i*2;
       A[i].z = B[i].z + (float) i*3;
    }
/* Copy back the results */ 
#pragma acc update self(A[0:size])
    printData(A, size);
    deleteData(A);
    deleteData(B);
    exit(0);
}

float3 * allocData(size_t size) {
    float3 * tmp;
    tmp = (float3 *) malloc(size*sizeof(float3));
/* Create the array on device. 
   Order matters.  The host copy must be allocated before 
   creating the device copy   */
#pragma acc enter data create(tmp[0:size])
    return tmp;
}

int deleteData(float3 * A) {
/* Delete the host copy.
   Order matters.  The device copy must be deleted before
   the host copy is freed.  */
#pragma acc exit data delete(A)
    free(A);
}

int initData(float3 *A, size_t size, float val) {
    size_t i;
    for (i=0; i < size; ++i) {
	A[i].x = val;
	A[i].y = val;
	A[i].z = val;
    }
/* Update the device with the initial values */
#pragma acc update device(A[0:size])
}

int printData(float3 *A, size_t size) {
    size_t i;
    printf("Values:\n");
    for (i=0; i < 10; ++i) {   
	printf("A[%d]=%f\n",i,A[i].x);
    } 
    printf("....\n");
    for (i=size-10; i < size; ++i) {   
	printf("A[%d]=%f\n",i,A[i].x);
    } 
}


