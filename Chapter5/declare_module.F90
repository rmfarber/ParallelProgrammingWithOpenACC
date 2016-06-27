#ifndef N
#define N 32 
#endif
#ifndef M
#define M 32 
#endif

module declare_data 

   integer, parameter :: dp = selected_real_kind(15, 307)
   real(kind=dp), dimension(:,:), allocatable :: A, B
!$acc declare create(A,B)

end module declare_data


module use_data

   use declare_data

contains

   subroutine fillData (j,val) 
!$acc routine vector
      integer(8), value :: j
      real(kind=dp), value :: val    
      integer(8) :: i
!$acc loop vector    
       do i=1,M
          A(j,i) = B(j,i) + ((j-1)*M)+(i-1) 
       end do

   end subroutine fillData
end module use_data

program declare_example
    use use_data

    integer(8) :: N1, M1
    integer(8) :: i,j
    N1=N
    M1=M

    allocate(A(N1,M1),B(N1,M1))
!$acc kernels
    B=2.5_dp
!$acc end kernels

!$acc parallel loop gang 
    do j=1,N1
        call fillData(j,2.5_dp)
    end do
!$acc update self(A)
    do j=2,N1,2
        print *, j, A(j,1), A(j,M1)
    enddo

    deallocate(A,B)

end program declare_example

