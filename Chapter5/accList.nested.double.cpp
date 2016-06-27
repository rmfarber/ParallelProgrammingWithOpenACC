#include <iostream>
#include "accList.h" 

using namespace std;
#ifndef N
#define N 16 
#endif

int main() {

    accList< accList<double> > A(N), B(N);
    for (int i=0; i < B.size(); ++i) {
	accList<double>*tmpA = new accList<double>(N);
	accList<double>*tmpB = new accList<double>(N);
	A.insert(i,tmpA);
	B.insert(i,tmpB);
        for (int j=0; j <= N; ++j) {
           B[i][j]=2.5;
        }
	delete tmpA;
	delete tmpB;
    } 
    B.accUpdateDevice();
    #pragma acc parallel loop collapse(2) present(A,B)
    for (int i=0; i < N; ++i) {
        for (int j=0; j < N; ++j) {
	   A[i][j]=B[i][j]+i;
        }
    } 
    A.accUpdateSelf();

    for(int i=0; i<5; ++i) {
	cout << "A[" << i << "]: " << A[i][0] << endl;
    }
    cout << "......" << endl;
    for(int i=N-5; i<N; ++i) {
	cout << "A[" << i << "]: " << A[i][0] << endl;
    }
}
    


