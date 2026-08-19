#!/bin/bash

#SBATCH --job-name="QE bands for ZrFCl"
#SBATCH --partition=short
#SBATCH --ntasks=32

ulimit -l  unlimited

module load materials/qe/7.2-openmpi

# Go to the band directory
cd $SLURM_SUBMIT_DIR

# Step 1: SCF — generates the converged charge density
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < ZrFCl.scf.in > ZrFCl.scf.out

# Step 2: Copy SCF charge density to dos/ BEFORE band NSCF overwrites out/
rm -rf ../dos/out
cp -r out ../dos/

# Step 3: Band NSCF — reads charge density from out/ along k-path
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < ZrFCl.nscf.in > ZrFCl.nscf.out

# Step 4: bands.x post-processing
bands.x < ZrFCl.band.in > ZrFCl.band.out
