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

#include <mpi.h>


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

    int rank = 0;
    int size = 1;

    //Initialize MPI and determine rank and size
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
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
    
#if _OPENACC
    acc_device_t device_type = acc_get_device_type();
    if ( acc_device_nvidia == device_type )
    {
        int ngpus=acc_get_num_devices(acc_device_nvidia);
        
        int devicenum=rank%ngpus;
        acc_set_device_num(devicenum,acc_device_nvidia);
    }
    // Call acc_init after acc_set_device_num to avoid multiple contexts on device 0 in multi GPU systems
    acc_init(device_type);
#endif /*_OPENACC*/
    #pragma acc enter data create(A,Aref,Anew,rhs)

    int ix_start = 1;
    int ix_end   = (NX - 1);

    // Ensure correctness if NY%size != 0
    int chunk_size = ceil( (1.0*NY)/size );

    int iy_start = rank * chunk_size;
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

    if ( rank == 0) printf("Jacobi relaxation Calculation: %d x %d mesh\n", NY, NX);

    if ( rank == 0) printf("Calculate reference solution and time serial execution.\n");
    StartTimer();
    poisson2d_serial( rank, iter_max, tol );
    double runtime_serial = GetTimer();
    
    //MPI Warm-up to establish CUDA IPC connections
    for (int i=0; i<2; ++i)
    {
        int top    = (rank == 0) ? (size-1) : rank-1;
        int bottom = (rank == (size-1)) ? 0 : rank+1;
        #pragma acc host_data use_device( A )
        {
            //1. Sent row iy_start (first modified row) to top receive lower boundary (iy_end) from bottom
            MPI_Sendrecv( &A[iy_start][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, top   , 0, &A[iy_end][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

            //2. Sent row (iy_end-1) (last modified row) to bottom receive upper boundary (iy_start-1) from top
            MPI_Sendrecv( &A[(iy_end-1)][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, bottom, 0, &A[(iy_start-1)][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, top   , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
    }

    //Wait for all processes to ensure correct timing of the parallel version
    MPI_Barrier( MPI_COMM_WORLD );
    if ( rank == 0) printf("Parallel execution.\n");
    StartTimer();
    int iter  = 0;
    real error = 1.0;
    
    #pragma acc update device(A[(iy_start-1):(iy_end-iy_start)+2][0:NX],rhs[iy_start:(iy_end-iy_start)][0:NX])
    while ( error > tol && iter < iter_max )
    {
        error = 0.0;

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
        
        #pragma acc kernels async(1)
        for( int ix = ix_start; ix < ix_end; ix++ )
        {
            A[iy_start][ix] = Anew[iy_start][ix];
            A[iy_end-1][ix] = Anew[iy_end-1][ix];
        }
        
        #pragma acc kernels async(2)
        for (int iy = iy_start+1; iy < iy_end-1; iy++)
        {
            for( int ix = ix_start; ix < ix_end; ix++ )
            {
                A[iy][ix] = Anew[iy][ix];
            }
        }
        
        real globalerror = 0.0;
        MPI_Allreduce( &error, &globalerror, 1, MPI_REAL_TYPE, MPI_MAX, MPI_COMM_WORLD );
        error = globalerror;
        
        #pragma acc wait(1)

        //Periodic boundary conditions
        int top    = (rank == 0) ? (size-1) : rank-1;
        int bottom = (rank == (size-1)) ? 0 : rank+1;
        #pragma acc host_data use_device( A )
        {
            //1. Sent row iy_start (first modified row) to top receive lower boundary (iy_end) from bottom
            MPI_Sendrecv( &A[iy_start][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, top   , 0, &A[iy_end][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

            //2. Sent row (iy_end-1) (last modified row) to bottom receive upper boundary (iy_start-1) from top
            MPI_Sendrecv( &A[(iy_end-1)][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, bottom, 0, &A[(iy_start-1)][ix_start], (ix_end-ix_start), MPI_REAL_TYPE, top   , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }
        
        #pragma acc kernels async(2)
        for (int iy = iy_start; iy < iy_end; iy++)
        {
                A[iy][0]      = A[iy][(NX-2)];
                A[iy][(NX-1)] = A[iy][1];
        }
        
        #pragma acc wait(2)
        
        if(rank == 0 && (iter % 100) == 0) printf("%5d, %0.6f\n", iter, error);
        
        iter++;
    }
    #pragma acc update self(A[(iy_start-1):(iy_end-iy_start)+2][0:NX])
    MPI_Barrier( MPI_COMM_WORLD );
    double runtime = GetTimer();

    if (check_results( rank, ix_start, ix_end, iy_start, iy_end, tol ) && rank == 0)
    {
        printf( "Num GPUs: %d.\n", size );
        printf( "%dx%d: 1 GPU: %8.4f s, %d GPUs: %8.4f s, speedup: %8.2f, efficiency: %8.2f%\n", NY,NX, runtime_serial/ 1000.0, size, runtime/ 1000.0, runtime_serial/runtime, runtime_serial/(size*runtime)*100 );
    }

    #pragma acc exit data delete(A,Aref,Anew,rhs)
    MPI_Finalize();
    return 0;
}

#include "poisson2d_serial.h"
