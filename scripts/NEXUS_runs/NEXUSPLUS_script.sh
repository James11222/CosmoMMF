#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J NEXUSPLUS_RUN
#SBATCH --mail-user=jamessunseri@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 35:00:00

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
cd /global/homes/j/james12/CosmoMMF/scripts/NEXUS_runs/
srun -n 1 -c 64 --cpu_bind=cores julia run_NEXUS_final.jl NEXUSPLUS
