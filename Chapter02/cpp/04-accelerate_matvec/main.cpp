/*
 *  Copyright 2014 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
#include <cstdlib>
#include <cstdio>
#include <omp.h>

#include "vector.h"
#include "vector_functions.h"
#include "matrix.h"
#include "matrix_functions.h"

#define N 200
#define MAX_ITERS 100
#define TOL 1e-12
int main() {
  vector x,b;
  vector r,p,Ap;
  matrix A;
  
  double one=1.0, zero=0.0;
  double normr, rtrans, oldtrans, p_ap_dot , alpha, beta;
  int iter=0;

  //create matrix
  allocate_3d_poisson_matrix(A,N);
    
  printf("Rows: %d, nnz: %d\n", A.num_rows, A.row_offsets[A.num_rows]);

  allocate_vector(x,A.num_rows);
  allocate_vector(Ap,A.num_rows);
  allocate_vector(r,A.num_rows);
  allocate_vector(p,A.num_rows);
  allocate_vector(b,A.num_rows);

  initialize_vector(x,100000);
  initialize_vector(b,1);
 

  waxpby(one, x, zero, x, p);
  matvec(A,p,Ap);
  waxpby(one, b, -one, Ap, r);
  
  rtrans=dot(r,r);
  normr=sqrt(rtrans);
  
  double st = omp_get_wtime();
  do {
    if(iter==0) {
      waxpby(one,r,zero,r,p);
    } else {
      oldtrans=rtrans;
      rtrans = dot(r,r);
      beta = rtrans/oldtrans;
      waxpby(one,r,beta,p,p);
    }
    
    normr=sqrt(rtrans);
  
    matvec(A,p,Ap);
    p_ap_dot = dot(Ap,p);

    alpha = rtrans/p_ap_dot;

    waxpby(one,x,alpha,p,x);
    waxpby(one,r,-alpha,Ap,r);

    if(iter%10==0)
      printf("Iteration: %d, Tolerance: %.4e\n", iter, normr);
    iter++;
  } while(iter<MAX_ITERS && normr>TOL);
  double et = omp_get_wtime();

  printf("Total Iterations: %d Total Time: %lfs\n", iter, (et-st));
  printf("Final Tolerance: %.4e\n", normr);

  free_vector(x);
  free_vector(r);
  free_vector(p);
  free_vector(Ap);
  free_matrix(A);

  return 0;
}
