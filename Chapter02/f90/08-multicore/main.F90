!
!  Copyright 2014 NVIDIA Corporation
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
program main
  use vector_mod
  use matrix_mod
  use omp_lib ! Used for timing
  implicit none
  real(8), allocatable :: x(:), b(:), r(:), p(:), Ap(:)
  type(matrix)         :: a
  real(8), parameter   :: one=1.0,       &
                          zero=0.0,      &
                          tol=1e-12
  integer, parameter   :: n=200,         &
                          max_iters=100 
  real(8)              :: normr, rtrans, oldtrans, p_ap_dot, alpha, beta
  integer              :: iter
  real                 :: st, et

  call allocate_3d_poisson_matrix(a,n)

  print *, "Rows:",a%num_rows,"nnz:",a%nnz

  call allocate_vector(x,a%num_rows)
  call allocate_vector(ap,a%num_rows)
  call allocate_vector(r,a%num_rows)
  call allocate_vector(p,a%num_rows)
  call allocate_vector(b,a%num_rows)

  call initialize_vector(x,100000d0)
  call initialize_vector(b,1d0)

  call waxpby(one, x, zero, x, p)
  call matvec(a,p,ap)
  call waxpby(one, b, -one, ap, r)

  rtrans=dot(r,r)
  normr=sqrt(rtrans)

  iter = 0
  st = omp_get_wtime()
  do while((iter.lt.max_iters).and.(normr.gt.tol)) 
    if(iter.eq.0) then
      call waxpby(one,r,zero,r,p)
    else
      oldtrans = rtrans
      rtrans = dot(r,r)
      beta = rtrans/oldtrans
      call waxpby(one,r,beta,p,p)
    endif

    normr = sqrt(rtrans)

    call matvec(a,p,ap)
    p_ap_dot = dot(ap,p)

    alpha = rtrans/p_ap_dot

    call waxpby(one,x,alpha,p,x)
    call waxpby(one,r,-alpha,ap,r)
    
    if (mod(iter,10).eq.0) then
      write(*,'(a,i3,a,es13.6)') "Iteration:",iter," Tolerance:",normr
    endif
    iter = iter + 1
  enddo
  et = omp_get_wtime()

  print *,"Total Iterations:",iter,"Time (s):",(et-st)
  
  call free_vector(x)
  call free_vector(ap)
  call free_vector(r)
  call free_vector(p)
  call free_vector(b)
  call free_matrix(a)
end program main
