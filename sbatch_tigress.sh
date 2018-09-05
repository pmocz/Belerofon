#!/bin/sh
#SBATCH -J B64
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH -t 0:00:10
#SBATCH -o output/OUTPUT.lsf
#SBATCH -e output/ERROR.lsf
#SBATCH --mail-user=philip.mocz@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=8000     # 80000 for 1024^3 run

module purge
module load intel
module load openmpi
module load fftw

./belerofon 40


