#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --mail-user=jamessunseri@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 8:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
module load python
conda activate baryon_env
srun -n 1 -c 64 --cpu_bind=cores python Generate_Halo_Hist.py ultralight