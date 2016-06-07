!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! module containing the output routine
MODULE m_io
 
  ! module dependencies
  USE m_config,  ONLY: nout, nx, ny, nz
  USE m_fields,  ONLY: qv

  IMPLICIT NONE

  CONTAINS

  !----------------------------------------------------------------------------  
  ! routine which outputs results (here only exemplary by computing a mean)
  SUBROUTINE write_output(ntstep)

    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN) :: ntstep ! current timestep

    ! local variables
    INTEGER :: i, j, k            ! loop indices
    REAL*8  :: qv_mean            ! scalar for computing the mean
    
    ! skip if this is not a output timestep
    IF (MOD(ntstep, nout) /= 0) RETURN
    
    !$acc data present(qv)

    ! compute mean of variable qv
    qv_mean = 0.0D0
    !$acc parallel 
    !$acc loop gang vector collapse(3) reduction(+:qv_mean)
    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          qv_mean = qv_mean + qv(i,j,k)
        END DO
      END DO
    END DO
    !$acc end parallel
    !$acc end data
    qv_mean = qv_mean / REAL(nx * ny * nz, KIND(qv_mean))

    ! echo result
    WRITE(*,"(A,I6,A,ES18.8)") "Step: ", ntstep, ", mean(qv) =", qv_mean

  END SUBROUTINE write_output

END MODULE m_io
