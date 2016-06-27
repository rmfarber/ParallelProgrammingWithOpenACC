#include <stdio.h>
#include <stdlib.h>

#ifndef N
#define N 32 
#endif
#ifndef M
#define M 32 
#endif

extern double A[N][M];
#pragma acc declare create (A)

#pragma acc routine vector
void setVal(int,double);
void printData();

