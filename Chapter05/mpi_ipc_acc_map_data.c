
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <openacc.h>

// CUDA runtime includes
#include <cuda_runtime_api.h>

#ifndef N
#define N 16
#endif

int main(int argc, char **argv) 
{
    int  rank, numprocs;
    MPI_Status status;
    double * devSharedMemory;
    double * sharedMemory;
    cudaIpcMemHandle_t mem_handle;

   
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    acc_set_device_num(0,acc_get_device_type());
    if (numprocs != 2) {
        printf("This example is written for exactly two ranks.\n");
	printf("Please run with -np 2. Aborting.\n");
	exit(1);
    }
    sharedMemory = (double*) malloc(N*sizeof(double));
    devSharedMemory = NULL;

    if ( rank == 0 ) {
	/* Rank 0 creates and gets an IPC handle to the device memory */
	#pragma acc enter data create(sharedMemory[0:N])
        cudaIpcGetMemHandle(&mem_handle,acc_deviceptr(sharedMemory));
    }	
	 
    if ( rank == 0 ) {
	/* Rank 0 sends the IPC Handle to the Rank 1 */
	MPI_Send(&mem_handle,CUDA_IPC_HANDLE_SIZE,MPI_BYTE,1,1,MPI_COMM_WORLD);	
    } else {
        /* Rank 1 gets the IPC Handle from the Rank 0 */
	MPI_Recv(&mem_handle,CUDA_IPC_HANDLE_SIZE,MPI_BYTE,0,1,MPI_COMM_WORLD,&status);	

	/* Open the IPC Handle and associate the memory to a local device pointer */
	cudaIpcOpenMemHandle(&devSharedMemory, mem_handle, cudaIpcMemLazyEnablePeerAccess);

	/* Map the Rank 1's shareMemory host array to the shared device memory. */	
	acc_map_data(sharedMemory,devSharedMemory,N*sizeof(double));
    }

    if ( rank == 1 ) {
	/* The Rank 1 sets the shared device memory values. */
	#pragma acc parallel loop present(sharedMemory) 
	for (size_t i=0; i < N; ++i) {
	   sharedMemory[i] = (double) i / (double) N;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if ( rank == 0 ) {
        /* The Rank 0 updates the host array and prints out the values set by the Rank 1. */
        #pragma acc update self(sharedMemory[0:N])
	for (size_t i=0; i < N; ++i) {
	    printf("sharedMemory[%d] = %f \n",i,sharedMemory[i]);
	} 
    }	
    MPI_Barrier(MPI_COMM_WORLD);

 
    if ( rank == 0 ) {
        /* Rank 0 deletes the device array */ 
	#pragma acc exit data delete(sharedMemory)
    } else {
	/* Rank 1 unmaps the data and closes the IPC Handle */
        acc_unmap_data(sharedMemory);
	cudaIpcCloseMemHandle(devSharedMemory);
    }	
    free(sharedMemory);
    MPI_Finalize();
    return 0;
}

            
