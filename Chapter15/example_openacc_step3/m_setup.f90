!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! module containing init and clean up routines
MODULE m_setup

  ! module dependencies
  USE m_config,  ONLY: nstop, nout, nx, ny, nz
  USE m_fields,  ONLY: t,qv
  USE m_timing,  ONLY: init_timers, start_timer, end_timer
  USE m_physics, ONLY: init_physics, finalize_physics
  IMPLICIT NONE

CONTAINS

  !----------------------------------------------------------------------------
  ! allocate memory and initialize modules and fields
  SUBROUTINE initialize()

    IMPLICIT NONE

    ! local variables
    INTEGER, PARAMETER :: itiminit = 1  ! timer ID
    INTEGER :: i, j, k                  ! loop indices

    ! print info
#ifdef _OPENACC
    WRITE(*,"(A)") "Running with OpenACC"
#else
    WRITE(*,"(A)") "Running without OpenACC"
#endif

    WRITE(*,"(A)") "Initialize"

    CALL init_timers()

    CALL start_timer( itiminit, "Initialization" )

    ! allocate memory
    ALLOCATE( t(nx,ny,nz), qv(nx,ny,nz) )

    ! allocate on the GPU
    !$acc enter data create(t,qv)

    ! initialize global fields
    DO k =1, nz
      DO j = 1, ny
        DO i = 1, nx
          t(i,j,k)  = 293.0D0 * (1.2D0 + 0.07D0 * COS(6.2D0 * REAL(i+j+k) / REAL(nx+ny+nz)))
          qv(i,j,k) = 1.0D-6 * (1.1D0 + 0.13D0 * COS(5.3D0 * REAL(i+j+k) / REAL(nx*ny*nz)))
        END DO
      END DO
    END DO

    ! initialize fields on the GPU
    !$acc update device(t,qv)

#ifdef _OPENACC
    ! call small kernel to initialize the GPU
    CALL initialize_gpu()
#endif
    ! initialize physics fields
    CALL init_physics()

    CALL end_timer( itiminit )

  END SUBROUTINE initialize

  !----------------------------------------------------------------------------
  SUBROUTINE initialize_gpu()

    IMPLICIT NONE

    ! local variables
    INTEGER :: temp(16)
    INTEGER :: i

    ! run a small kernel before any timings to make sure
    ! the GPU is initialized

    !$acc parallel loop
    DO i = 1, 16
      temp(i) = 1
    END DO

    IF (SUM(temp) == 16) THEN
      WRITE(*,"(A)") "GPU initialized"
    ELSE
      WRITE(*,"(A,I4)") "Error: Problem encountered initializing the GPU"
      STOP
    END IF

  END SUBROUTINE initialize_gpu

  !----------------------------------------------------------------------------
  ! deallocate memory
  SUBROUTINE cleanup()

    IMPLICIT NONE
    
    ! deallocate on the GPU
    !$acc exit data delete(t,qv)

    DEALLOCATE( t, qv )

    CALL finalize_physics()

  END SUBROUTINE cleanup

END MODULE m_setup
