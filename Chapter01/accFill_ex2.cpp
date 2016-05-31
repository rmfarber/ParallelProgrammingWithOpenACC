#include <iostream>
using namespace std;
#include <cstdlib>
#include <cassert>

int main(int argc, char *argv[])
{
  if(argc < 2) {
    cerr << "Use: nCount" << endl;
    return -1;
  }

  int nCount = atoi(argv[1])*1000000;

  if(nCount < 0) {
    cerr << "ERROR: nCount must be greater than zero!" << endl;
    return -1;
  }

   // allocate variables
   char *status = new char[nCount];

   // Here is where we fill the status vector
   int sum=0; // Important to zero-initialize sum
#pragma acc kernels create( status[0:nCount] ) 
   {
     for(int i=0; i < nCount; i++)
       status[i] = 1;
#pragma acc loop vector reduction(+:sum)
     for(int i=0; i < nCount; i++)
       sum += status[i];
   }
   
   cout << "Final sum is " << (sum/1000000) << " millions" << endl;

   // Sanity check to see that array is filled with ones
   assert(sum == nCount);
   
   delete [] status;
   
   return 0; // normal exit
}
