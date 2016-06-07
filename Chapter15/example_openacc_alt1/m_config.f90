!This file is released under terms of MIT license`
!See LICENSE.txt for more information

! module containing run configuration
MODULE m_config

  ! module variables
  INTEGER :: nx    = 128  ! domain size in the longitudinal direction 
  INTEGER :: ny    = 128  ! domain size in the latitudinal direction 
  INTEGER :: nz    = 60   ! domain size in the vertical direction
  INTEGER :: nstop = 100  ! number of time steps
  INTEGER :: nout  = 20   ! output ever nout steps

END MODULE m_config
