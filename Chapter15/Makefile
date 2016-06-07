# Makefile for OpenACC example code

# Default
COMPILER=cray
NV=nvcc
# Adapt architecture to GPU
NVFLAGS=-arch=sm_35

# Set compiler and options
ifeq ($(COMPILER),pgi)
 FC=pgf90
 LD=$(FC)
 ifeq ($(TARGET),gpu) 
  FFLAGS=-O3 -Mpreprocess -ta=nvidia
  LFLAGS=-ta=nvidia
  ifdef WITHCUDA
   LFLAGS+=-Mcuda
  endif
 else
  FFLAGS=-O3 -Mpreprocess -mp
  LFLAGS=-mp
 endif
endif

ifeq ($(COMPILER),cray)
 FC=ftn
 LD=$(FC)
 ifeq ($(TARGET),gpu) 
  FFLAGS=-O3 -eZ -hacc
  LFLAGS=-hacc
 else
  FFLAGS=-O3 -eZ -homp -hnoacc
  LFLAGS=-homp
 endif
endif


ifeq ($(COMPILER),gnu)
 FC=gfortran
 LD=$(FC)
 FFLAGS=-O3 -cpp -fopenmp
 LFLAGS=-fopenmp
endif


# Object files
OBJ=main.o m_config.o m_fields.o m_io.o m_parameterizations.o m_physics.o m_setup.o m_timing.o

ifdef WITHCUDA
 CUDAOBJ=saturation_adjustment.o
endif

# Rules
.PHONY: example_serial example_openmp example_openacc_step1 example_openacc_step2 example_openacc_step3 example_openacc_alt1 example_openacc_alt2 clean

# compile in example_name/ directory
all: example_serial example_openmp example_openacc_step1 example_openacc_step2 example_openacc_step3 example_openacc_alt1 example_openacc_alt2

example_serial:
	@make -C $@ -f ../Makefile exe TARGET=cpu EXE=$@

example_openmp:
	@make -C $@ -f ../Makefile exe TARGET=cpu EXE=$@

example_openacc_step1:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc_step2:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc_step3:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc_alt1:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc_alt2:
	@make -C $@ -f ../Makefile exe TARGET=gpu WITHCUDA=1 EXE=$@

exe : $(OBJ) $(CUDAOBJ)
	@$(LD) $(LFLAGS) $(OBJ) $(CUDAOBJ) -o $(EXE)

clean : 
	rm -f */*.cub */*.i */*.o */*.mod */*.ptx */*~

# Suffix Rules 
.SUFFIXES: .o .mod .f90

%.o : %.f90
	@echo "compiling $<"
	@$(FC) -c $(FFLAGS) -o $@ $<

saturation_adjustment.o : saturation_adjustment.cu
	@echo "compiling $<"
	$(NV) -c $(NVFLAGS) saturation_adjustment.cu

# dependencies
m_io.o: m_config.o m_fields.o
m_physics.o: m_config.o m_fields.o m_parameterizations.o
m_setup.o: m_config.o m_fields.o m_timing.o
main.o: m_config.o m_io.o m_physics.o m_setup.o m_timing.o
ifdef WITHCUDA
 m_parameterizations.o: saturation_adjustment.o
endif
