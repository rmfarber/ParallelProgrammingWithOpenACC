program mandelbrot_main
use mandelbrot_mod
use openacc
implicit none
integer      :: num_blocks
integer(1)   :: image(HEIGHT, WIDTH)
integer      :: iy, ix
integer      :: block, block_size, block_start
integer      :: starty, endy
real         :: startt, stopt
character(8) :: arg

call acc_init(acc_device_nvidia)

call getarg(1, arg)
num_blocks = 8
if ( arg /= '' ) then
  read(arg, '(I10)') num_blocks
endif
print *,'num_blocks',num_blocks
block_size = (HEIGHT*WIDTH)/num_blocks

image = 0

call cpu_time(startt)
!$acc data create(image(HEIGHT,WIDTH))
do block=0,(num_blocks-1)
  starty = block  * (WIDTH/num_blocks) + 1
  endy   = min(starty + (WIDTH/num_blocks), WIDTH)
  !$acc parallel loop async(mod(block,2))
  do iy=starty,endy
    do ix=1,HEIGHT
      image(ix,iy) = min(max(int(mandelbrot(ix-1,iy-1)),0),MAXCOLORS)
    enddo
  enddo
  !$acc update self(image(:,starty:endy)) async(mod(block,2))
enddo
!$acc wait
!$acc end data
call cpu_time(stopt)

print *,"Time:",(stopt-startt)

call write_pgm(image,'image.pgm')
end
