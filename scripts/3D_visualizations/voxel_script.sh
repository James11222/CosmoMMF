#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J Voxel_RUN
#SBATCH --mail-user=jamessunseri@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 5:00:00

#OpenMP settings:
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
cd /global/homes/j/james12/CosmoMMF/scripts/
module load python
conda activate Baryon_Env
srun -n 1 -c 64 --cpu_bind=cores python compute_voxel.py
