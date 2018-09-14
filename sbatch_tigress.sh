#!/bin/sh
#SBATCH -J bfn512
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH -t 24:00:00
#SBATCH -o output/OUTPUT.lsf
#SBATCH -e output/ERROR.lsf
#SBATCH --mail-user=philip.mocz@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=24000     # 80000 for 1024^3 run

export OMP_NUM_THREADS=40

module purge
module load intel
module load openmpi
module load fftw

srun -u -n 1  --tasks-per-node=40 ./belerofon 40


