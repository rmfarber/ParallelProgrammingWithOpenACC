!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! module containing the driving code of the 
! physical parameterizations

MODULE m_physics

  ! module dependencies
  USE m_config, ONLY: nx, ny, nz
  USE m_fields, ONLY: t, qv
  USE m_parameterizations, ONLY: saturation_adjustment, microphysics

  IMPLICIT NONE

  REAL*8, ALLOCATABLE :: qc(:,:,:)

CONTAINS

  !----------------------------------------------------------------------------
  ! driving routine for the physical parameterizations
  SUBROUTINE physics()

    IMPLICIT NONE
 
    ! call a first physical parameterization
    CALL saturation_adjustment(nx, ny, nz, t, qc, qv)

    ! call a second physical parameterization  
    CALL microphysics(nx, ny, nz, t, qc, qv)

  END SUBROUTINE physics

  !-----------------------------------------------------------------------------
  ! Initialization of local physics array
  SUBROUTINE init_physics()
   
    IMPLICIT NONE

    ! allocate memory
    ALLOCATE( qc(nx,ny,nz) )

    !$acc enter data create(qc)

  END SUBROUTINE init_physics

  !----------------------------------------------------------------------------
  ! Finalize the local physics array
  SUBROUTINE finalize_physics()
    
    IMPLICIT NONE
    
    ! deallocate memory
    !$acc exit data delete(qc)

    DEALLOCATE(qc)

   END SUBROUTINE finalize_physics

END MODULE m_physics
