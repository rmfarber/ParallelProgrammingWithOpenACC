#!/bin/bash
#SBATCH --job-name=openacc_example
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --gres=gpu:1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=job.out

# Export env variables
export OMP_NUM_THREADS=1
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=536870912
 
# Set this to avoid segmentation faults
ulimit -s unlimited
ulimit -a

echo "====================================================="
echo "============== JOB OUTPUT BEGINS ===================="
echo "====================================================="

export G2G=1
export MV2_USE_CUDA=1
export MV2_USE_GPUDIRECT=0
export COSMO_NPROC_NODEVICE=0 #no async io currently

srun -u -n 1 example_openacc1/example_openacc1

echo "====================================================="
echo "============== JOB OUTPUT ENDS ===================="
echo "====================================================="

