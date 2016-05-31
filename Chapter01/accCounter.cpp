#include <iostream>
using namespace std;
#include <cstdlib>
#include <cassert>
#include <chrono>
using namespace std::chrono;

int main(int argc, char *argv[])
{
   if(argc < 2) {
      cerr << "Use: nCount" << endl;
      return 1;
      }
    int nCount = atoi(argv[1]);

   if(nCount <= 0) {
      cerr << "ERROR: nCount must be greater than zero!" << endl;
      return 1;
   }

   high_resolution_clock::time_point t1 = high_resolution_clock::now();

   // Here is where we define and increment the counter.
    int counter=0;
#pragma acc parallel loop 
   for( int i=0; i < nCount; i++) {
#ifdef USE_ATOMIC
     #pragma acc atomic update
#endif
     counter++;
   }

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

   cout << "counter " << counter 
	<< " Duration " << time_span.count() << " second" << endl;
   assert(counter == nCount);

   return 0; // normal exit
}
