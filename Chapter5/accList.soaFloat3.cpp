
#include <iostream>
#include "soaFloat3.h"
#include "accList.h"
#ifndef N
#define N 6
#define M 1024
#endif

int main() {

  accList<soaFloat3> C{N};
  for (size_t j = 0; j < C.size(); ++j) {
    soaFloat3 tmpC{M};
    C.insert(j, tmpC);
  }

  accList<soaFloat3> A{N}, B{N};
  for (size_t j = 0; j < A.size(); ++j) {
    soaFloat3 tmpA{M};
    soaFloat3 tmpB{M};
    A.insert(j, tmpA);
    B.insert(j, tmpB);
    B[j].setValue(2.5, 3.5, 4.5);
  }
  A.accUpdateDevice();
  B.accUpdateDevice();
  #pragma acc parallel loop present(A,B)
  for (size_t j = 0; j < N; ++j) {
    A[j].setValue(2.5, 3.5, 4.5);
    A[j].add(B[j]);
  }
  A.accUpdateSelf();
  for (size_t j = 0; j < 5; ++j) {
    cout << "A[" << j << "]: " << A[j].x[1] << "," 
	 << A[j].y[1] << "," << A[j].z[1] << endl;
  }
  cout << "......" << endl;
  for (size_t j = N-5; j < N; ++j) {
    cout << "A[" << j << "]: " << A[j].x[1] << "," 
	 << A[j].y[1] << "," << A[j].z[1] << endl;
  }
}
