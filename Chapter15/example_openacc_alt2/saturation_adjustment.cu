#include  <stdio.h>
#include "cuda.h"

extern "C" {

	void saturation_adjustment_cuda(int ntot,
                                    double *t, double *qc, double *qv,
                                    double cs1, double cs2, double cs3, double cs4, double t0);

}

//--------------------------------------
// saturation adjustment kernel
__global__  void saturation_adjustment_kernel(int ntot,
                                              double *t, double *qc, double *qv,
                                              double cs1, double cs2, double cs3, double cs4, double t0)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if ( tid < ntot )
    {
        qv[tid] = qv[tid] + cs1 * exp( cs2 * ( t[tid] - t0 ) / ( t[tid] - cs3 ) );
        qc[tid] = cs4 * qv[tid];
    }
}

//--------------------------------------
// CUDA routine calling the saturation adjustment kernel
 void saturation_adjustment_cuda(int ntot,
                                 double *t, double *qc, double *qv,
                                 double cs1, double cs2, double cs3, double cs4, double t0)
{
    // set CUDA grid dimensions
    const int THREADS_PER_BLOCK = 128; // number of gpu threads per block
    const int NUMBER_OF_BLOCKS = ceil(ntot / THREADS_PER_BLOCK); 

    // calling CUDA kernel
    saturation_adjustment_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(ntot, t, qc, qv,
                                                                          cs1, cs2, cs3, cs4, t0);
}

