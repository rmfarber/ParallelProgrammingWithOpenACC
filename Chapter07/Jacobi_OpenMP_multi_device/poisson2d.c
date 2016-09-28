/* Copyright (c) 2016, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENACC
#include <openacc.h>
#endif /*_OPENACC*/
#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/

#include "common.h"

#define NY 4096
#define NX 4096

real A[NY][NX];
real Aref[NY][NX];
real Anew[NY][NX];
real rhs[NY][NX];

int main(int argc, char** argv)
{
    int iter_max = 1000;
    
    const real tol = 1.0e-5;

    int thread_num = 0;
    int num_threads = 1;
    
    real globalerror = 0.0;
    
#pragma omp parallel firstprivate(thread_num,num_threads) shared(globalerror)
{
#ifdef _OPENMP
    thread_num = omp_get_thread_num();
    num_threads = omp_get_num_threads();
#endif /*_OPENMP*/
#pragma omp master
{
    // set rhs
    for (int iy = 1; iy < NY-1; iy++)
    {
        for( int ix = 1; ix < NX-1; ix++ )
        {
            const real x = -1.0 + (2.0*ix/(NX-1));
            const real y = -1.0 + (2.0*iy/(NY-1));
            rhs[iy][ix] = expr(-10.0*(x*x + y*y));
        }
    }
}
    
#if _OPENACC
    acc_device_t device_type = acc_get_device_type();
    if ( acc_device_nvidia == device_type )
    {
        int ngpus=acc_get_num_devices(acc_device_nvidia);
        
        int devicenum=thread_num%ngpus;
        acc_set_device_num(devicenum,acc_device_nvidia);
    }
    acc_init(device_type);
#endif /*_OPENACC*/

    int ix_start = 1;
    int ix_end   = (NX - 1);

    // Ensure correctness if NY%num_threads != 0
    int chunk_size = ceil( (1.0*NY)/num_threads );

    int iy_start = thread_num * chunk_size;
    int iy_end   = iy_start + chunk_size;

    // Do not process boundaries
    iy_start = max( iy_start, 1 );
    iy_end = min( iy_end, NY - 1 );
    
    //OpenACC Warm-up
    #pragma acc kernels
    for( int iy = 0; iy < NY; iy++)
    {
        for( int ix = 0; ix < NX; ix++ )
        {
            Aref[iy][ix] = 0.0;
            A[iy][ix] = 0.0;
        }
    }

#pragma omp master
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", NY, NX);

#pragma omp master
    printf("Calculate reference solution and time serial execution.\n");
    double runtime_serial = 0.0;
#pragma omp master
{
    StartTimer();
    poisson2d_serial( iter_max, tol );
    runtime_serial = GetTimer();
}
    //Wait for all threads to ensure correct timing of the parallel version
#pragma omp barrier
#pragma omp master
{
    printf("Parallel execution.\n");
    StartTimer();
}
    int iter  = 0;
    real error = 1.0;
    
    #pragma acc data copy(A[(iy_start-1):(iy_end-iy_start)+2][0:NX]) copyin(rhs[iy_start:(iy_end-iy_start)][0:NX]) create(Anew[iy_start:(iy_end-iy_start)][0:NX])
    while ( error > tol && iter < iter_max )
    {
        error = 0.0; 
#pragma omp single
        globalerror = 0.0;
#pragma omp barrier

        #pragma acc kernels
        for (int iy = iy_start; iy < iy_end; iy++)
        {
            for( int ix = ix_start; ix < ix_end; ix++ )
            {
                Anew[iy][ix] = -0.25 * (rhs[iy][ix] - ( A[iy][ix+1] + A[iy][ix-1]
                                                       + A[iy-1][ix] + A[iy+1][ix] ));
                error = fmaxr( error, fabsr(Anew[iy][ix]-A[iy][ix]));
            }
        }
        
#pragma omp critical
        globalerror = fmaxr(globalerror,error);
#pragma omp barrier
        error = globalerror;
        
        #pragma acc kernels
        for (int iy = iy_start; iy < iy_end; iy++)
        {
            for( int ix = ix_start; ix < ix_end; ix++ )
            {
                A[iy][ix] = Anew[iy][ix];
            }
        }

        //Periodic boundary conditions and halo updates
        #pragma acc update self( A[iy_start:1][0:NX], A[(iy_end-1):1][0:NX] )
#pragma omp barrier
        if ( 0 == (iy_start-1) )
        {
            for( int ix = 1; ix < NX-1; ix++ )
            {
                A[0][ix]      = A[(NY-2)][ix];
            }
        }
        if ( NY-1 == iy_end )
        {        
            for( int ix = 1; ix < NX-1; ix++ )
            {
                    A[(NY-1)][ix] = A[1][ix];
            }
        }
        #pragma acc update device( A[(iy_start-1):1][0:NX], A[iy_end:1][0:NX] )

        #pragma acc kernels
        for (int iy = iy_start; iy < iy_end; iy++)
        {
                A[iy][0]      = A[iy][(NX-2)];
                A[iy][(NX-1)] = A[iy][1];
        }
        
#pragma omp master
{
        if((iter % 100) == 0) printf("%5d, %0.6f\n", iter, error);
}
        
        iter++;
    }

#pragma omp barrier
    double runtime = 0.0;
#pragma omp master
    runtime = GetTimer();

    if (check_results( thread_num, ix_start, ix_end, iy_start, iy_end, tol ))
    {
#pragma omp master
{
        printf( "Num GPUs: %d.\n", num_threads );
        printf( "%dx%d: 1 GPU: %8.4f s, %d GPUs: %8.4f s, speedup: %8.2f, efficiency: %8.2f%\n", NY,NX, runtime_serial/ 1000.0, num_threads, runtime/ 1000.0, runtime_serial/runtime, runtime_serial/(num_threads*runtime)*100 );
}
    }

} //#pragma omp parallel
    return 0;
}

#include "poisson2d_serial.h"
