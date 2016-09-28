program mandelbrot_main
use mandelbrot_mod
implicit none
!integer, parameter   :: NUM_BLOCKS=8
integer(1) :: image(HEIGHT, WIDTH)
integer :: iy, ix
real :: startt, stopt
image = 0

call cpu_time(startt)
!$acc parallel loop
do iy=1,width
  do ix=1,HEIGHT
    image(ix,iy) = min(max(int(mandelbrot(ix-1,iy-1)),0),MAXCOLORS)
  enddo
enddo
call cpu_time(stopt)

print *,"Time:",(stopt-startt)

call write_pgm(image,'image.pgm')
end
