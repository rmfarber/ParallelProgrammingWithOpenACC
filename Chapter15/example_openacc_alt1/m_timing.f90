!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! simple timing utilities to time code sections
MODULE m_timing

  IMPLICIT NONE

  ! module constants
  INTEGER, PARAMETER :: ntimer=10            ! number of timers

  ! module variables
  REAL*8             :: rtimer(ntimer)       ! timers
  CHARACTER(32)      :: timertag(ntimer)     ! unique timer tags
  INTEGER            :: icountold(ntimer), & ! tick (start of timer section)
                        icountrate,        & ! countrate of SYSTEM_CLOCK()
                        icountmax            ! maximum counter value of SYSTEM_CLOCK()

  CONTAINS

  !----------------------------------------------------------------------------  
  ! initialize timers
  SUBROUTINE init_timers()

    IMPLICIT NONE

    ! initialize all timers
    rtimer(:)   = 0.0D0
    timertag(:) = ""
    icountold(:) = 0

    ! initialize count rate and max
    CALL SYSTEM_CLOCK( COUNT_RATE=icountrate, COUNT_MAX=icountmax )

  END SUBROUTINE init_timers

  !----------------------------------------------------------------------------
  ! start timer for this id
  SUBROUTINE start_timer(id, tag)

    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN)       :: id
    CHARACTER(*), INTENT(IN) :: tag

    ! check if timer ID is within legal range
    IF (id < 1 .OR. id > ntimer) THEN
      WRITE(*,"(A,I4,A,I4)") "Error: timer id=", id, "exceeds maximum timer number", ntimer
      STOP
    END IF

    ! check that timer has not alreay been started
    IF (LEN_TRIM(timertag(id)) /= 0) THEN
      WRITE(*,"(A,I4)") "Error: timer already started previously, id:", id
      STOP
    END IF

    ! check that tag is non-empty
    IF (LEN_TRIM(tag) == 0) THEN
      WRITE(*,"(A,I4)") "Error: empty tag provided, id:", id
      STOP
    END IF

    ! save tag
    timertag(id) = TRIM(tag)

    ! start timer
    !$acc wait
    CALL SYSTEM_CLOCK( COUNT=icountold(id) )
    
  END SUBROUTINE start_timer

  !----------------------------------------------------------------------------
  ! end timer compute elapsed time
  SUBROUTINE end_timer(id)

    IMPLICIT NONE

    ! arguments
    INTEGER, INTENT(IN) :: id

    ! local variables
    INTEGER             :: icountnew

    ! check if timer ID is within legal range
    IF (id < 1 .OR. id > ntimer) THEN
      WRITE(*,"(A,I4,A,I4)") "Error: timer id=", id, "exceed max timer number", ntimer
      STOP
    END IF

    ! check if timer has been started previously
    IF (LEN_TRIM(timertag(id)) == 0) THEN
      WRITE(*,"(A,I4)") "Error: Need to call start_timer before end_timing, id:", id
      STOP
    END IF

   ! get current time
   !$acc wait
   CALL SYSTEM_CLOCK( COUNT=icountnew )

   ! compute elapsed time
   rtimer(id) = ( REAL(icountnew - icountold(id), KIND(rtimer(id))) ) /   &
                REAL(icountrate, KIND(rtimer(id)))

  END SUBROUTINE end_timer

  !----------------------------------------------------------------------------
  ! echo timers to stdout
  SUBROUTINE print_timers()

    IMPLICIT NONE

    ! local variables
    INTEGER :: id

    ! print all non-zero timers in microsecond
    WRITE(*,"(A)") "----------------------------"
    WRITE(*,"(A)") "Timers:"
    WRITE(*,"(A)") "----------------------------"
    DO id = 1, ntimer
      IF ( rtimer(id) > 0.0D0 ) THEN
        WRITE(*,"(A15,A2,F8.2,A)") timertag(id), ": ", rtimer(id)*1.0D3, " ms"
      END IF
    END DO    
    WRITE(*,"(A)") "----------------------------"

  END SUBROUTINE print_timers
  
END MODULE m_timing
