#include "declare_global_static.h"

double A[N][M];

int main() {
   int i;
#pragma acc parallel loop gang 
   for(i=0;i<N;++i) {
     setVal(i,2.5);
   }
#pragma acc update self(A[0:N][0:M])
   printData();
   exit(0);
}


