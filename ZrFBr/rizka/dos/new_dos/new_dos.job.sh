#!/bin/bash

#SBATCH --job-name="ZrFBr dos+pdos"
#SBATCH --partition=short
#SBATCH --ntasks=32
#SBATCH --output=slurm_ZrFBr_new_dos%j.log
#SBATCH --error=slurm_ZrFBr_new_dos%j.err

# NOTE: Runs SCF then dense NSCF, then computes total DOS and projected DOS (PDOS).

ulimit -l unlimited

module load materials/qe/7.2-openmpi

# Go to the dos directory
cd $SLURM_SUBMIT_DIR

# Step 1: SCF to generate charge density
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x <ZrFBr.scf.in >ZrFBr.scf.out

sleep 5

# Step 2: Dense NSCF for DOS
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x <ZrFBr.nscf.in >ZrFBr.nscf.out

sleep 5

# Step 3: Compute total DOS
dos.x <ZrFBr.dos.in >ZrFBr.dos.out

sleep 5

# Step 4: Compute projected DOS (PDOS)
projwfc.x <ZrFBr.pdos.in >ZrFBr.pdos.out
