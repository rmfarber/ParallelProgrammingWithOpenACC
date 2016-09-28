COMPILER ?= PGI-tesla
CC = mpicc
CFLAGS = -DUSE_DOUBLE
ifeq ($(COMPILER),GCC)
	CFLAGS += -march=native -O3 -lm
else ifeq ($(COMPILER),PGI-tesla)
	CFLAGS += -Minfo=accel -fast -acc -ta=tesla
else ifeq ($(COMPILER),PGI-multicore)
	CFLAGS += -Minfo=accel -fast -acc -ta=multicore
endif
all: run

poisson2d: poisson2d.c common.h poisson2d_serial.h Makefile
	$(CC) $(CFLAGS) poisson2d.c -o poisson2d

clean:
	rm -f poisson2d poisson2d.*.nvprof 

run: poisson2d
	mpirun -np 4 ./poisson2d

profile: poisson2d
	mpirun -np 4 nvprof --process-name "MPI %q{OMPI_COMM_WORLD_RANK}" --context-name "MPI %q{OMPI_COMM_WORLD_RANK}" -o poisson2d.%q{OMPI_COMM_WORLD_RANK}.nvprof ./poisson2d

