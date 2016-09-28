program mandelbrot_main
use mandelbrot_mod
use omp_lib
use openacc
implicit none
integer, parameter   :: NUM_BLOCKS=64
integer(1) :: image(HEIGHT, WIDTH)
integer :: iy, ix
integer :: block, block_size, block_start
integer :: starty, endy
integer :: num_gpus, my_gpu, queue
real :: startt, stopt

num_gpus = acc_get_num_devices(acc_device_nvidia)

! This region is used to absorb the start-up cost
! so that it doesn't skew timing results
!$omp parallel num_threads(num_gpus)
call acc_init(acc_device_nvidia)
call acc_set_device(omp_get_thread_num(),acc_device_nvidia)
!$omp end parallel

block_size = (HEIGHT*WIDTH)/NUM_BLOCKS

image = 0
queue = 1

!$omp parallel num_threads(num_gpus) private(my_gpu) firstprivate(queue)
my_gpu = omp_get_thread_num()
call acc_set_device_num(my_gpu,acc_device_nvidia)
print *, "Thread:",my_gpu,"is using GPU",acc_get_device_num(acc_device_nvidia)

startt = omp_get_wtime()
!$acc data create(image(HEIGHT,WIDTH))
!$omp do schedule(static,1)
do block=0,(num_blocks-1)
  starty = block  * (WIDTH/NUM_BLOCKS) + 1
  endy   = min(starty + (WIDTH/NUM_BLOCKS), WIDTH)
  !$acc parallel loop async(queue)
  do iy=starty,endy
    do ix=1,HEIGHT
      image(ix,iy) = min(max(int(mandelbrot(ix-1,iy-1)),0),MAXCOLORS)
    enddo
  enddo
  !$acc update self(image(:,starty:endy)) async(queue)
  queue = mod((queue+1),2)
enddo
!$acc wait
!$acc end data
!$omp end parallel
stopt = omp_get_wtime()

print *,"Time:",(stopt-startt)

call write_pgm(image,'image.pgm')
end
