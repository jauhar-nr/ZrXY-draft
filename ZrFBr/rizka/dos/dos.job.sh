#!/bin/bash

#SBATCH --job-name="QE dos for ZrFBr"
#SBATCH --partition=short
#SBATCH --ntasks=32

# NOTE: Make sure to run the band job first to generate the charge density 'out' directory.

ulimit -l  unlimited

module load materials/qe/7.2-openmpi

# Go to the dos directory
cd $SLURM_SUBMIT_DIR

# Run NSCF for DOS
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < ZrFBr.nscf.in > ZrFBr.nscf.out
dos.x < ZrFBr.dos.in > ZrFBr.dos.out
