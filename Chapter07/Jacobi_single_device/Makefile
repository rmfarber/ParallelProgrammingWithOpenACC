COMPILER ?= PGI-tesla
CFLAGS = -DUSE_DOUBLE
ifeq ($(COMPILER),GCC)
	CC = gcc
	CFLAGS += -march=native -O3 -lm
else ifeq ($(COMPILER),PGI-tesla)
	CC = pgcc
	CFLAGS += -Minfo=accel -fast -acc -ta=tesla
else ifeq ($(COMPILER),PGI-multicore)
	CC = pgcc
	CFLAGS += -Minfo=accel -fast -acc -ta=multicore
endif
all: run

poisson2d: poisson2d.c common.h Makefile
	$(CC) $(CFLAGS) poisson2d.c -o poisson2d

clean:
	rm -f poisson2d poisson2d.*.nvprof 

run: poisson2d
	./poisson2d

profile: poisson2d
	nvprof -o poisson2d.nvprof ./poisson2d
