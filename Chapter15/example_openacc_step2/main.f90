!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! Main code. This represent a simplified atmospheric model 
! structure with only the physical parameterization
! and the output.

PROGRAM main

  ! module dependencies
  USE m_config,  ONLY: nstop
  USE m_physics, ONLY: physics
  USE m_io,      ONLY: write_output
  USE m_setup,   ONLY: initialize, cleanup
  USE m_timing,  ONLY: start_timer, end_timer, print_timers
  
  IMPLICIT NONE

  ! local variables
  INTEGER :: ntstep
  INTEGER, parameter :: itimloop = 5
  
  CALL initialize()
  
  !----------------------------------------------------------------------------
  ! time loop

  WRITE(*,"(A)") "Start of time loop"

  CALL start_timer(itimloop, "Time loop")

  DO ntstep = 1, nstop

    ! call the physical parameterizations
    CALL physics()
    
    ! call outputs
    CALL write_output( ntstep )

  END DO
  
  CALL end_timer( itimloop )

  WRITE(*,"(A)") "End of time loop"

  !----------------------------------------------------------------------------

  CALL print_timers()

  CALL cleanup()

END PROGRAM main
