#include <iostream>
#include <cstdlib>
#include <cstdint>
#include "accList.h" 
using namespace std;
#ifndef N
#define N 1024
#endif

int main() {

    accList<double> A(N), B(N);
    for (int i=0; i < B.size(); ++i) {
        B[i]=2.5;
    } 
    B.accUpdateDevice();
    #pragma acc parallel loop gang vector present(A,B)
    for (int i=0; i < A.size(); ++i) {
	   A[i]=B[i]+i;
    } 
    A.accUpdateSelf();

    for(int i=0; i<10; ++i) {
	cout << "A[" << i << "]: " << A[i] << endl;
    }
    cout << "......" << endl;
    for(int i=N-10; i<N; ++i) {
	cout << "A[" << i << "]: " << A[i] << endl;
    }
}
    


