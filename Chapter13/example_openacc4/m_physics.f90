!This file is released under terms of BSD license`
!See LICENSE.txt for more information

! module containing the driving code of the 
! physical parameterizations

MODULE m_physics

  ! module dependencies
  USE m_config, ONLY: nx, ny, nz
  USE m_fields, ONLY: t, qv 
  USE m_parametrizations, ONLY: saturation_adjustment, microphysics

  IMPLICIT NONE

  REAL*8, ALLOCATABLE :: qc(:,:,:)

CONTAINS

  !----------------------------------------------------------------------------
  ! driving routine for the physical parameterizations
  SUBROUTINE physics()

    IMPLICIT NONE

    INTEGER :: i,j    ! loop indices
    !$acc data present(t,qc,qv) 
    !$acc parallel
    !$acc loop gang
    DO j = 1,ny
      !$acc loop vector
      DO i = 1,nx
        ! call a first physical parametrization
        CALL saturation_adjustment(nz, i, j, t(:,:,:), qc(:,:,:), qv(:,:,:))

        ! call a second physical parametrization  
        CALL microphysics(nz, i, j, t(:,:,:), qc(:,:,:), qv(:,:,:))
      END DO
   END DO
   !$acc end parallel
   !$acc end data
  END SUBROUTINE physics
  !----------------------------------------------------------------------------
  ! Initialization of local physics array
  SUBROUTINE init_physics()

    IMPLICIT NONE

    ! allocate memory
    ALLOCATE ( qc(nx,ny,nz) )

    !$acc enter data create(qc)
    !
 END SUBROUTINE init_physics

  !----------------------------------------------------------------------------
  ! Finalize the local physics array
  SUBROUTINE finalize_physics()

    IMPLICIT NONE

    ! deallocate memory
    !$acc exit data delete(qc)

    DEALLOCATE (qc)

 END SUBROUTINE finalize_physics

END MODULE m_physics
