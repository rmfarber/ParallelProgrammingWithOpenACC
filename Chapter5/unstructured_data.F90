#ifndef N
#define N 32 
#endif
#ifndef M
#define M 32 
#endif

module unstruct

   integer, parameter :: dp = selected_real_kind(15, 307)
   real(kind=dp), dimension(:,:), allocatable :: A, B

contains

   subroutine allocateData (N1,M1)
       integer(8) :: N1,M1 
       allocate(A(N1,M1),B(N1,M1))
       !$acc enter data create(A,B)
   end subroutine allocateData

   subroutine deallocateData()
       !$acc exit data delete(A,B)
       deallocate(A,B)
   end subroutine deallocateData

end module unstruct

program unstruct_example
    use unstruct

    integer(8) :: N1, M1
    integer(8) :: i,j
    N1=N
    M1=M

    call allocateData(N1,M1)
    B=2.5_dp
!$acc update device(B)
!$acc kernels present(A,B)
    do j=1,N1
       do i=1,M1
          A(j,i) = B(j,i) + ((j-1)*M1)+(i-1) 
       end do
    end do
!$acc end kernels   
!$acc update self(A)
    do j=2,N1,2
        print *, j, A(j,1), A(j,M1)
    enddo

    call deallocateData()

end program unstruct_example

