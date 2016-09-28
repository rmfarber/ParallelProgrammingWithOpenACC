# Makefile for OpenACC example code

# Default
COMPILER=cray

# Set compiler and options
ifeq ($(COMPILER),pgi)
 FC=pgf90
 LD=$(FC)
 ifeq ($(TARGET),gpu) 
  FFLAGS=-O3 -Mpreprocess -ta=nvidia
  LFLAGS=-ta=nvidia
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
OBJ=main.o m_config.o m_fields.o m_io.o m_parametrizations.o m_physics.o m_setup.o m_timing.o

# Rules
.PHONY: example_serial example_openmp example_openacc1 example_openacc2 example_openacc3 example_openacc4 clean

# compile in example_name/ directory
example_serial:
	@make -C $@ -f ../Makefile exe TARGET=cpu EXE=$@

example_openmp:
	@make -C $@ -f ../Makefile exe TARGET=cpu EXE=$@

example_openacc1:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc2:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc3:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

example_openacc4:
	@make -C $@ -f ../Makefile exe TARGET=gpu EXE=$@

exe : $(OBJ)
	@$(LD) $(LFLAGS) $(OBJ) -o $(EXE)

clean : 
	rm -f */*.cub */*.i */*.o */*.mod */*.ptx */*~

# Suffix Rules 
.SUFFIXES: .o .mod .f90

%.o : %.f90
	@echo "compiling $<"
	@$(FC) -c $(FFLAGS) -o $@ $<

# dependencies
m_io.o: m_config.o m_fields.o
m_physics.o: m_config.o m_fields.o m_parametrizations.o
m_setup.o: m_config.o m_fields.o m_timing.o
main.o: m_config.o m_io.o m_physics.o m_setup.o m_timing.o

