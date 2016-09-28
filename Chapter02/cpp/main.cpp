#include <cstdlib>
#include <cstdio>

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
  allocate_3d_poission_matrix(A,N);
    
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

  printf("Total Iterations: %d\n", iter);

  free_vector(x);
  free_vector(r);
  free_vector(p);
  free_vector(Ap);
  free_matrix(A);

  return 0;
}
