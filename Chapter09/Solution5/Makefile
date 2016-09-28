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
APP := 3Dwave_acc_multicore

UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
   NX := 256
   NY := 256
   NZ := 256
   NSTEPS := 200
   DEFINES := WIN32
   OBJ := obj
   EXE := exe
   CCFLAGS := -D$(DEFINES) -fast -Minline -Minfo=all -acc -ta=multicore -tp=haswell
   INCLUDES := -I. -I/c/"Program Files"/PGI/win64/2015/cuda/7.5/include
   LIBS := -L/c/"Program Files"/"NVIDIA GPU Computing Toolkit"/CUDA/v7.5/lib/x64 cudart.lib
   NUMA :=
else
   NX := 944
   NY := 944
   NZ := 944
   NSTEPS := 500
   CCFLAGS := -D$(DEFINES) -fast -Minline -Minfo=all -acc -ta=multicore -tp=haswell
   LIBS := -L$(CUDA_ROOT)/lib64 -lcudart
   NUMA := numactl --cpunodebind=0 --physcpubind=0,1,2,3,4,5,6,7
endif

all: build run verify

build: $(APP).c
	$(CC) $(INCLUDES) $(CCFLAGS) -o $(APP).$(EXE) $(APP).c $(LIBS)

run: ./$(APP).$(EXE)
	$(NUMA) ./$(APP).$(EXE) $(NX) $(NY) $(NZ) $(NSTEPS)

verify:	../Verify/verify.$(EXE)
	../Verify/verify.$(EXE) $(NX) $(NSTEPS) wave3d.xline.ref wave3d.xline

clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.$(OBJ) *.dwf *.pdb prof
