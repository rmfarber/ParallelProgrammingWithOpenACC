! Example demonstrating cuda ipc communication in fortran with OpenACC
! This example shows how one per-GPU master rank can create a shared
! memory buffer on the GPU and send a handle to a slave which can use
! this handle to directly work on the GPU memory without memcopies.
! Obviously this is meant for ranks that use the same physical GPU.
! Author C.Angerer

module cuda_ipc_api
    use iso_c_binding
    use openacc
    implicit none

    !data types and size taken from cuda driver_types.h
    integer, parameter :: CUDA_IPC_HANDLE_SIZE = 64
    integer, parameter :: cudaIpcMemLazyEnablePeerAccess = 1
    type, bind(c) :: cudaIpcMemHandle_t
        character(kind=c_char) :: opaque(CUDA_IPC_HANDLE_SIZE)
    end type cudaIpcMemHandle_t

    interface    
    !cudaError_t cudaIpcGetMemHandle(cudaIpcMemHandle_t *handle, void *devPtr);
    integer(c_int) function cudaIpcGetMemHandle(handle, dev_ptr) bind(c,name="cudaIpcGetMemHandle")
        import
        !ignore type kind and rank for dev_ptr; still requires "device"
        !attribute
        !dir$ ignore_tkr (tkr) dev_ptr 
        type(cudaIpcMemHandle_t) :: handle
        type(c_devptr), device :: dev_ptr(*)
    end function cudaIpcGetMemHandle

    !cudaError_t cudaIpcOpenMemHandle(void **devPtr, cudaIpcMemHandle_t handle, unsigned int flags);
    !flags must be cudaIpcMemLazyEnablePeerAccess
    integer(c_int) function cudaIpcOpenMemHandle(dev_ptr, handle, flags) bind(c,name="cudaIpcOpenMemHandle")
        import
        type(c_devptr), intent(out) :: dev_ptr
        type(cudaIpcMemHandle_t), value :: handle
        integer(c_int), value :: flags
    end function cudaIpcOpenMemHandle

    !cudaError_t cudaIpcCloseMemHandle (void *devPtr)
    integer(c_int) function cudaIpcCloseMemHandle(dev_ptr) bind(c,name="cudaIpcCloseMemHandle")
        import
        type(c_devptr), intent(in), value :: dev_ptr
    end function cudaIpcCloseMemHandle
    end interface
end module cuda_ipc_api

program mpi_test_gpu
    use mpi
    use cuda_ipc_api
    use openacc
    use iso_c_binding
    integer, allocatable :: sharedMemory(:)
    integer:: N, ierr, rank, num_procs, status(MPI_Status_size)
    type(cudaIpcMemHandle_t) :: mem_handle
    type(c_devptr) :: slave_dev_mem_ptr, host_dev_mem_ptr

    !Note: use Mpi_comm_split_type with MPI_COMM_TYPE_SHARED to create a
    !node-local communicator. This works on Titan because every node only has
    !one GPU. For multi-GPU nodes use some app-specific way to determine 
    !what ranks actually run on the same GPU.
    call MPI_Init (ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    if(num_procs .ne. 2) then
        print *, "This example is written for exactly two ranks."
        print *, "Please run with -np 2. Aborting."
        stop
    endif
    N=4

    !Note: we are allocating space for the whole shared data in slave host memory;
    !Only the device memory is mapped into one
    !luckily this makes portability easier when compiling without CUDA IPC: just add
    !mpi send/receives and use the available host buffers
    allocate (sharedMemory(N))
    sharedMemory(:) = -1

    !Step 1: Master rank prepares shared memory and creates a handle
    if ( rank == 0) then
        print *, "Master creating shared memory handle"
        !$acc enter data create(sharedMemory)
        !$acc host_data use_device(sharedMemory)
        ierr = cudaIpcGetMemHandle(mem_handle, sharedMemory)
        !$acc end host_data
    end if

    !Step 2: exchange handle to shared memory between master and slave and
    !associate received device pointer with OpenACC variable
    if ( rank == 0) then
        print *, "Master sending shared memory handle"
        call MPI_Send(mem_handle,CUDA_IPC_HANDLE_SIZE,MPI_CHAR,1,0,MPI_COMM_WORLD, ierr)
    else
        call MPI_Recv(mem_handle,CUDA_IPC_HANDLE_SIZE,MPI_CHAR,0,0,MPI_COMM_WORLD,status, ierr)
        print *, "Slave received shared memory handle; mapping it to host buffer"
        ierr = cudaIpcOpenMemHandle(slave_dev_mem_ptr, mem_handle, cudaIpcMemLazyEnablePeerAccess)
        call acc_map_data(sharedMemory, slave_dev_mem_ptr, sizeof(sharedMemory))
    end if

    !Step 3: work on the shared GPU memory in the slave
    if (rank == 1) then
        print *, "Slave working on shared GPU memory"
        !$acc kernels
        sharedMemory(:) = 42
        !$acc end kernels
    end if

    !Step 4: synchronize
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !Step 5: on master: check that we see the changes that slave did in sharedMemory
    if ( rank == 0) then
        print *, "Master reading shared GPU memory (should all be 42)"
        !$acc update self(sharedMemory)
        print *, sharedMemory
    end if

    !Step 6: synchronize
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    !Step 7: clean up
    if ( rank == 0) then
        print *, "Master cleaning up"
        !acc exit data create(sharedMemory)
    else
        print *, "Slave cleaning up"
        call acc_unmap_data(sharedMemory)
        ierr = cudaIpcCloseMemHandle(slave_dev_mem_ptr)
    end if
    deallocate (sharedMemory)

    call MPI_Finalize ( ierr )
end program mpi_test_gpu
