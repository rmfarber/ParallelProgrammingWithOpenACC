#!/bin/bash

module load pgi/16.4

for i in 2 3 4 5; do
   echo Solution$i
   cd Solution$i; make build; cd ..
done
