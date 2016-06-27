
CC=pgcc
FC=pgfortran
CXX=pgc++ 
NVCC=nvcc
MPICC=mpicc
MPIF90=mpif90
MPIRUN=mpirun
OPT=-fast 
ACC=-acc  -Minfo=accel
CXXOPT=--c++11 -w -I.
CUDA=-Mcuda
NVOPT=-O3

EXES= accList_soaFloat3 accList_double accList_nested_double \
      acc_malloc acc_map array_of_structs declare_module \
      declare_global_static jagged_array myList unstructured_data2D \
      unstructured_data_c unstructured_data_f unstructured_data_struct \
      mpi_ipc_acc_map_data_c mpi_ipc_acc_map_data_f

all: $(EXES)

acc_malloc: acc_malloc.out
	-./acc_malloc.out

acc_malloc.out: acc_malloc.c cudaSet.cu
	$(NVCC) $(NVOPT) -c -o cudaSet.o cudaSet.cu
	$(CC) $(OPT) $(CUDA) $(ACC) -o acc_malloc.out cudaSet.o acc_malloc.c

acc_map: acc_map.out
	-./acc_map.out

acc_map.out: acc_map.c
	$(CC) $(OPT) $(ACC) -o acc_map.out acc_map.c
		
array_of_structs: array_of_structs.out
	-./array_of_structs.out

array_of_structs.out: array_of_structs.c
	$(CC) $(OPT) $(ACC) -o array_of_structs.out array_of_structs.c
	
declare_global_static: declare_global_static.out
	-./declare_global_static.out
	
declare_global_static.out: declare_global_static.c declare_global_static_routine.c 
	$(CC) $(OPT) $(ACC) -o declare_global_static.out declare_global_static.c declare_global_static_routine.c
		
declare_module: declare_module.out
	-./declare_module.out

declare_module.out: declare_module.F90
	$(FC) $(OPT) $(ACC) -o declare_module.out declare_module.F90

accList_soaFloat3: accList.soaFloat3.out
	-./accList.soaFloat3.out

accList.soaFloat3.out: accList.soaFloat3.cpp soaFloat3.h accList.h
	$(CXX) $(CXXOPT) $(OPT) $(ACC) -o accList.soaFloat3.out accList.soaFloat3.cpp

accList_double: accList.double.out
	-./accList.double.out

accList.double.out: accList.double.cpp accList.h
	$(CXX) $(CXXOPT) $(OPT) $(ACC) -o accList.double.out accList.double.cpp

accList_nested_double: accList.nested.double.out
	-./accList.nested.double.out

accList.nested.double.out: accList.nested.double.cpp accList.h
	$(CXX) $(CXXOPT) $(OPT) $(ACC) -o accList.nested.double.out accList.nested.double.cpp

jagged_array: jagged_array.out
	-./jagged_array.out

jagged_array.out: jagged_array.c
	$(CC) $(OPT) $(ACC) -o jagged_array.out jagged_array.c

mpi_ipc_acc_map_data_c: mpi_ipc_acc_map_data_c.out
	-$(MPIRUN) -np 2 ./mpi_ipc_acc_map_data_c.out

mpi_ipc_acc_map_data_c.out: mpi_ipc_acc_map_data.c
	$(MPICC) $(OPT) $(ACC) $(CUDA) -o mpi_ipc_acc_map_data_c.out mpi_ipc_acc_map_data.c

mpi_ipc_acc_map_data_f: mpi_ipc_acc_map_data_f.out
	-$(MPIRUN) -np 2 ./mpi_ipc_acc_map_data_f.out

mpi_ipc_acc_map_data_f.out: mpi_ipc_acc_map_data.f90
	$(MPIF90) $(OPT) $(ACC) $(CUDA) -o mpi_ipc_acc_map_data_f.out mpi_ipc_acc_map_data.f90

myList: myList.out
	-./myList.out

myList.out: myList.cpp myList.h
	$(CXX) $(CXXOPT) $(OPT) $(ACC) -o myList.out myList.cpp

unstructured_data2D: unstructured_data2D.out
	-./unstructured_data2D.out

unstructured_data2D.out: unstructured_data2D.c
	$(CC) $(OPT) $(ACC) -o unstructured_data2D.out unstructured_data2D.c

unstructured_data_struct: unstructured_data_struct.out
	-./unstructured_data_struct.out

unstructured_data_struct.out: unstructured_data.struct.c
	$(CC) $(OPT) $(ACC) -o unstructured_data_struct.out unstructured_data.struct.c

unstructured_data_c: unstructured_data_c.out
	-./unstructured_data_c.out

unstructured_data_c.out: unstructured_data.c
	$(CC) $(OPT) $(ACC) -o unstructured_data_c.out unstructured_data.c

unstructured_data_f: unstructured_data_f.out
	-./unstructured_data_f.out

unstructured_data_f.out: unstructured_data.F90
	$(FC) $(OPT) $(ACC) -o unstructured_data_f.out unstructured_data.F90

clean:
	rm -f *.out *.o *.mod
