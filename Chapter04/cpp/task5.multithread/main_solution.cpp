#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include <openacc.h>
#include "mandelbrot.h"
#include "constants.h"

using namespace std;

int main() {
  
  size_t bytes=WIDTH*HEIGHT*sizeof(unsigned int);
  unsigned char *image=(unsigned char*)malloc(bytes);
  int num_blocks=64, block_size = (HEIGHT/num_blocks)*WIDTH;
  FILE *fp=fopen("image.pgm","wb");
  fprintf(fp,"P5\n%s\n%d %d\n%d\n","#comment",WIDTH,HEIGHT,MAX_COLOR);

  int num_gpus = acc_get_num_devices(acc_device_nvidia);
// This parallel section eats the cost of initializing the devices to
// prevent the initialization time from skewing the results.
#pragma omp parallel num_threads(num_gpus)
{
  acc_init(acc_device_nvidia);
  acc_set_device_num(omp_get_thread_num(),acc_device_nvidia);
}
  printf("Found %d NVIDIA GPUs.\n", num_gpus);

  double st = omp_get_wtime();
#pragma omp parallel num_threads(num_gpus)
{
  int queue = 1;
  int my_gpu = omp_get_thread_num();
  acc_set_device_num(my_gpu,acc_device_nvidia);
  printf("Thread %d is using GPU %d\n", my_gpu, acc_get_device_num(acc_device_nvidia));
#pragma acc data create(image[WIDTH*HEIGHT])
{
  #pragma omp for schedule(static,1)
  for(int block = 0; block < num_blocks; block++ ) {
    int start = block * (HEIGHT/num_blocks),
        end   = start + (HEIGHT/num_blocks);
#pragma acc parallel loop async(queue)
    for(int y=start;y<end;y++) {
      for(int x=0;x<WIDTH;x++) {
        image[y*WIDTH+x]=mandelbrot(x,y);
      }
    }
#pragma acc update self(image[block*block_size:block_size]) async(queue)
    queue = (queue + 1) % 2; 
  }
}
#pragma acc wait
} // OMP Parallel 
  
  double et = omp_get_wtime();
  printf("Time: %lf seconds.\n", (et-st));
  fwrite(image,sizeof(unsigned char),WIDTH*HEIGHT,fp);
  fclose(fp);
  free(image);
  return 0;
}
