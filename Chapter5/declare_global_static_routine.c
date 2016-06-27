/*************************************************************
OpenACC Book
Example of use "declare create" with a global static array
**************************************************************/

#include "declare_global_static.h"
#pragma acc routine vector
void setVal(int i, double val) {
    int j;
#pragma acc loop vector
    for(j=0;j<M;++j) {
	A[i][j] = val + (j*M)+i;
    }
}

void printData() {
    size_t i,j,m;
    m=M-1;
    printf("Values:\n");
    for (i=0; i < 5; ++i) {
        printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i][0],i,m,A[i][m]);
    }
    printf("....\n");
    for (i=N-5; i < N; ++i) {
        printf("A[%d][0]=%f A[%d][%d]=%f\n",i,A[i][0],i,m,A[i][m]);
    }
}



