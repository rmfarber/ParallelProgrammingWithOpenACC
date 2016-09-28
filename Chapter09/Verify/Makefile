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

UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
   DEFINES := WIN32
   OBJ := obj
   EXE := exe
   CCFLAGS := -D$(DEFINES) -fast -Minline
   LIBS := 
else
   CCFLAGS := -D$(DEFINES) -fast -Minline
   LIBS := 
endif


all: build run

build: verify.c
	$(CC) $(INCLUDES) $(CCFLAGS) -o verify.$(EXE) verify.c $(LIBS)

run: ./verify.$(EXE)
	./run.sh $^

clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.$(OBJ) *.dwf *.pdb prof
