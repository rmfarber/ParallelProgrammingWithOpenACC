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

#ifndef LAPLACE2D_SERIAL_H
#define LAPLACE2D_SERIAL_H

int global_result_correct = 1;

void poisson2d_serial( int iter_max, real tol )
{
    int iter  = 0;
    real error = 1.0;
    #pragma acc data copy(Aref) copyin(rhs) create(Anew)
    while ( error > tol && iter < iter_max )
    {
        error = 0.0;

#pragma acc kernels
        for( int iy = 1; iy < NY-1; iy++)
        {
            for( int ix = 1; ix < NX-1; ix++ )
            {
                Anew[iy][ix] = -0.25 * (rhs[iy][ix] - ( Aref[iy][ix+1] + Aref[iy][ix-1]
                                                       + Aref[iy-1][ix] + Aref[iy+1][ix] ));
                error = fmaxr( error, fabsr(Anew[iy][ix]-Aref[iy][ix]));
            }
        }

#pragma acc kernels
        for( int iy = 1; iy < NY-1; iy++)
        {
            for( int ix = 1; ix < NX-1; ix++ )
            {
                Aref[iy][ix] = Anew[iy][ix];
            }
        }

        //Periodic boundary conditions
#pragma acc kernels
        for( int ix = 1; ix < NX-1; ix++ )
        {
                Aref[0][ix]     = Aref[(NY-2)][ix];
                Aref[(NY-1)][ix] = Aref[1][ix];
        }
#pragma acc kernels
        for( int iy = 1; iy < NY-1; iy++ )
        {
                Aref[iy][0]     = Aref[iy][(NX-2)];
                Aref[iy][(NX-1)] = Aref[iy][1];
        }

        if((iter % 100) == 0) printf("%5d, %0.6f\n", iter, error);

        iter++;
    }
}

int check_results( int thread_num, int ix_start, int ix_end,  int iy_start, int iy_end, real tol )
{
    int result_correct = 1;
    for( int iy = iy_start; iy < iy_end && (result_correct == 1); iy++)
    {
        for( int ix = ix_start; ix < ix_end && (result_correct == 1); ix++ )
        {
            if ( fabs ( Aref[iy][ix] - A[iy][ix] ) >= tol )
            {
                fprintf(stderr,"[Thread%d] ERROR: A[%d][%d] = %f does not match %f (reference)\n", thread_num, iy,ix, A[iy][ix], Aref[iy][ix]);
                result_correct = 0;
            }
        }
    }
#pragma omp critical
{
    global_result_correct = min(global_result_correct,result_correct);
}
#pragma omp barrier
    result_correct = global_result_correct;
    return result_correct;
}

#endif // LAPLACE2D_SERIAL_H
