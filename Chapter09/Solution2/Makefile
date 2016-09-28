################################################################################
#
# Copyright 2014 NVIDIA Corporation.  All rights reserved.
#
# Please refer to the NVIDIA end user license agreement (EULA) associated
# with this source code for terms and conditions that govern your use of
# this software. Any use, reproduction, disclosure, or distribution of
# this software and related documentation outside the terms of the EULA
# is strictly prohibited.
#
################################################################################

CC := pgcc
DEFINES := FP32
INCLUDES := -I.
OBJ := o
EXE := out
APP := 3Dwave_omp_dynamic

UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
   NX := 256
   NY := 256
   NZ := 256
   NSTEPS := 500
   DEFINES := WIN32
   OBJ := obj
   EXE := exe
else
   NX := 944
   NY := 944
   NZ := 944
   NSTEPS := 500
endif

CCFLAGS := -D$(DEFINES) -fast -Minline -Minfo=all,ccff -mp

all: build run verify

build: $(APP).c
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(APP).$(EXE) $(APP).c

run:
	make run8
	make run4
	make run1

run8: ./$(APP).$(EXE)
	export OMP_NUM_THREADS=8; ./$(APP).$(EXE) $(NX) $(NY) $(NZ) $(NSTEPS)

run4: ./$(APP).$(EXE)
	export OMP_NUM_THREADS=4; ./$(APP).$(EXE) $(NX) $(NY) $(NZ) $(NSTEPS)

run1: ./$(APP).$(EXE)
	export OMP_NUM_THREADS=1; ./$(APP).$(EXE) $(NX) $(NY) $(NZ) $(NSTEPS)
	cp wave3d.xline wave3d.xline.ref
	cp wave3d.xline.ref ../Solution3
	cp wave3d.xline.ref ../Solution4
	cp wave3d.xline.ref ../Solution5

verify:	../Verify/verify.$(EXE)
	../Verify/verify.$(EXE) $(NX) $(NSTEPS) wave3d.xline.ref wave3d.xline

clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.$(OBJ) *.dwf *.pdb prof
