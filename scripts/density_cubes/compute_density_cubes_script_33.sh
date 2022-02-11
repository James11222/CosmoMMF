#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --mail-user=jamessunseri@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 15:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
source /global/common/software/m3035/conda-activate.sh 3.7
module load h5py-parallel
cd /global/u2/j/james12/CosmoMMF/scripts/
srun -n 1 -c 64 --cpu_bind=cores python create_density_cubes.py 33 all
