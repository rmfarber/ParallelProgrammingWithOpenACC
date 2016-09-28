!
!  Copyright 2014 NVIDIA Corporation
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module vector_mod
  implicit none
  public :: initialize_vector,allocate_vector,free_vector
  public :: dot, waxpby
  contains
  subroutine initialize_vector(vector,value)
    implicit none
    real(8), intent(out) :: vector(:)
    real(8), intent(in)  :: value
    vector(:) = value
  end subroutine
  subroutine allocate_vector(vector,length)
    implicit none
    real(8), allocatable :: vector(:)
    integer              :: length
    allocate(vector(length))
  end subroutine allocate_vector
  subroutine free_vector(vector)
    implicit none
    real(8), allocatable :: vector(:)
    deallocate(vector)
  end subroutine
  function dot(x, y)
    implicit none
    real(8), intent(in) :: x(:), y(:)
    real(8)             :: dot, tmpsum
    integer             :: i, length

    length = size(x)
    tmpsum = 0.0
    do i=1,length
      tmpsum = tmpsum + x(i)*y(i)
    enddo

    dot = tmpsum
  end function
  subroutine waxpby(alpha, x, beta, y, w)
    implicit none
    real(8), intent(in)  :: alpha, beta, x(:), y(:)
    real(8), intent(out) :: w(:)
    integer             :: i, length

    length = size(x)
    !$acc kernels
    do i=1,length
      w(i) = alpha*x(i) + beta*y(i)
    enddo
    !$acc end kernels

  end subroutine
end module vector_mod
