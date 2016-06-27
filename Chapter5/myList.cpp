#ifndef N
#define N 1024
#endif
#include <iostream>
#include <string.h>
#include "myList.h" 

int main() {

    myList<double> A(N), B(N);
    for (int i=0; i < B.size(); ++i) {
        B[i]=2.5;
    } 
    B.accUpdateDevice();
    #pragma acc parallel loop present(A,B)
    for (int i=0; i < A.size(); ++i) {
	   A[i]=B[i]+i;
    } 
    A.accUpdateSelf();
    for(int i=0; i<10; ++i) {
	cout << "A[" << i << "]: " << A[i] << endl;
    }
    exit(0);
}
    


