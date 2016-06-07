!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! module containing the physical parameterizations
MODULE m_parameterizations

  USE iso_c_binding, ONLY: C_LOC

  IMPLICIT NONE

  ! module constants
  REAL*8 :: cs1 = 1.0D-6, cs2 = 0.02D0, cs3 = 7.2D0, cs4 = 0.1D0, t0 = 273.0D0
  REAL*8 :: cm1 = 1.0D-6, cm2 = 25.0D0, cm3 = 0.2D0, cm4 = 100.0D0

CONTAINS

  !----------------------------------------------------------------------------
  ! Physical parameterization example code 1
  ! Note: The code here do not reflect actual physical equations.
  !       The code only reflects typical computational pattern encountered in the
  !       physics of an atmospheric model
  SUBROUTINE saturation_adjustment(npx, npy, nlev, t, qc, qv)

    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN)    :: npx, npy, nlev  ! dimensions of inputs
    REAL*8,  INTENT(IN)    :: t(:,:,:)        ! temperature
    REAL*8,  INTENT(OUT)   :: qc(:,:,:)       ! cloud water content
    REAL*8,  INTENT(INOUT) :: qv(:,:,:)       ! water vapour content

    ! Interface to CUDA wrapper function using iso_c_binding
    INTERFACE
      ! Pass pointer and dimensions of ijk fields to dycore
      SUBROUTINE saturation_adjustment_cuda( ntot, t, qc, qv,             &
                                             cs1, cs2, cs3, cs4, t0     ) &
           BIND(c, name='saturation_adjustment_cuda')
        USE, INTRINSIC :: iso_c_binding
        INTEGER(C_INT), VALUE      :: ntot
        TYPE(C_PTR), VALUE         :: t, qc, qv
        REAL(KIND=C_DOUBLE), VALUE :: cs1, cs2, cs3, cs4, t0
      END SUBROUTINE saturation_adjustment_cuda
    END INTERFACE

    ! Call the CUDA wrapper routine 

    !$acc data present(t,qv,qc)

    !$acc host_data use_device(t,qv,qc)
    CALL saturation_adjustment_cuda(npx*npy*nlev,             &
                                    C_LOC(t(1,1,1)),          &
                                    C_LOC(qc(1,1,1)),         &
                                    C_LOC(qv(1,1,1)),         &
                                    cs1, cs2, cs3, cs4, t0 )
    !$acc end host_data

    !$acc end data

  END SUBROUTINE saturation_adjustment

  !----------------------------------------------------------------------------
  ! Physical parameterization example code 2
  ! Note: The code here do not reflect actual physical equations.
  !       The code only reflects typical computational pattern encountered in the
  !       physics of an atmospheric model
  SUBROUTINE microphysics(npx, npy, nlev, t, qc, qv)

    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN)   :: npx, npy, nlev  ! dimension of inputs
    REAL*8, INTENT(INOUT) :: t(:,:,:)        ! temperature
    REAL*8, INTENT(IN)    :: qc(:,:,:)       ! cloud water content
    REAL*8, INTENT(INOUT) :: qv(:,:,:)       ! water vapour content

    ! local variables
    INTEGER :: i, j, k   ! loop indices

    !$acc data present(t,qv,qc)

    ! do the computation
    !$acc parallel
    !$acc loop seq
    DO k = 2, nlev
      !$acc loop gang 
      DO j = 1, npy
        !$acc loop vector
        DO i = 1, npx
          qv(i, j, k) = qv(i,j,k-1) + cm1*(t(i,j,k)-cm2)**cm3
          t(i, j, k)  = t(i, j, k)*( 1.0D0 - cm4*qc(i,j,k)+qv(i,j,k) )
        END DO
      END DO
    END DO
    !$acc end parallel

    !$acc end data

  END SUBROUTINE microphysics

END MODULE m_parameterizations
