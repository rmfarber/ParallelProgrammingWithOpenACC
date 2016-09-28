
Example source from chapter on Advanced Data Management

Type "make" to build and run all examples.

Code tested using PGI compiler version 16.5.
The "acc_malloc" use both pgcc and nvcc.
The "mpi_ipc_acc_map_data" examples use mpicc or mpif90

To run an individual example: make <example>

accList_soaFloat3: C++ list class using a soaFloat3 class as data member
accList_double: C++ list class using a fundamental data type as data member
accList_nested_double: C++ list class using another list class
acc_malloc: Using acc_malloc to declare device data passed to a CUDA kernel 
acc_map: Create a device memory pool and mapping host data to this pool 
array_of_structs: Create an array of stucts on the device
declare_module: Using Fortran module data in a device routine
declare_global_static: Using a C global static array in a device routine
jagged_array: Creating a jagged array on the device 
myList: Simple C++ list class
unstructured_data2D: Using unstructured data region with 2 dimensional array 
unstructured_data_c: Using unstructured data region in C
unstructured_data_f: Using unstructured data region in Fortran 
unstructured_data_struct: Using unstructured data region with a struct
mpi_ipc_acc_map_data_c: Sharing device data across MPI processes. C version
mpi_ipc_acc_map_data_f: Sharing device data across MPI processes. Fortran version
