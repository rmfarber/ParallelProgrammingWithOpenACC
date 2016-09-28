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

CONTAINS

  !----------------------------------------------------------------------------
  ! driving routine for the physical parameterizations
  SUBROUTINE physics()

    IMPLICIT NONE
 
    ! local variables
    REAL*8 :: qc(nx,ny,nz)   ! temporary variable used in the physics
   
    ! call a first physical parametrization
    CALL saturation_adjustment(nx, ny, nz, t, qc, qv)

    ! call a second physical parametrization  
    CALL microphysics(nx, ny, nz, t, qc, qv)

  END SUBROUTINE physics

END MODULE m_physics
