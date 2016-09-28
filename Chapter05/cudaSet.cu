/****************************************
ADD Description
TODO
*****************************************/


#include <cuda.h>

__global__ void setVal(double * B, size_t size, double val)
{
    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(; tid < size; tid += stride)
        B[tid] = val;
}

extern "C" 
{
void cudaSet(double * B, size_t size, double val) {
    setVal<<<1,128>>>(B,size,val);
}
}
