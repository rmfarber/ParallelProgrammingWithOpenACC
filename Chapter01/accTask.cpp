#include <iostream>
using namespace std;
#include <cstdlib>
#include <cassert>
#include <chrono>
using namespace std::chrono;

// a sequential task 
// a task that runs in a thread
#pragma acc routine worker
inline char task(int index, int nLoop)
{
  long counter=0;
  for(long i=0; i < 1000000*nLoop; i++)
      counter += (i>=0)?1:0;

  // return 1 if the counter is correct
  return( ((counter/1000000) == nLoop)?1:0 );
}

int main(int argc, char *argv[])
{
  if(argc < 3) {
    cerr << "Use: nCount nLoop" << endl;
    return -1;
  }
  
  int nCount = atoi(argv[1]);
  int nLoop = atoi(argv[2]);
  
  if(nCount < 0 || nLoop < 0) {
    cerr << "ERROR: both nCount and nLoop must be greater than zero!" << endl;
    return -1;
  }
  
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  // Here is where we evaluate the task(s)
  int sum=0;
#pragma acc parallel loop reduction(+:sum)
  for(int i=0; i < nCount; i++)
    sum += task(i,nLoop);
  
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast< duration<double> >(t2 - t1);
  
  cout << "Duration " << time_span.count() << " second" << endl;
  cout << "Final sum is " << sum << endl;
  
  // Sanity check to see that array is filled with ones
  assert(sum == nCount);
  
  return 0; // normal exit
}
