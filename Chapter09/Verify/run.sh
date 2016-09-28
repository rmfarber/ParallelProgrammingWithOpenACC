#!/bin/bash

EXECNAME=$1

NX=256
NSTEPS=2000

echo +++ Compare Results +++
if [ -e wave3d.xline.ref ] && [ -e wave3d.xline ]; then
   ./$EXECNAME $NX $NSTEPS
fi
echo ""
