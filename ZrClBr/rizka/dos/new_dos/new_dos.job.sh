#!/bin/bash

#SBATCH --job-name="ZrClB dos+pdos"
#SBATCH --partition=short
#SBATCH --ntasks=32
#SBATCH --output=slurm_ZrClbr_new_dos%j.log
#SBATCH --error=slurm_ZrClbr_new_dos%j.err

# NOTE: Runs SCF then dense NSCF, then computes total DOS and projected DOS (PDOS).

ulimit -l unlimited

module load materials/qe/7.2-openmpi

# Go to the dos directory
cd $SLURM_SUBMIT_DIR

# Step 1: SCF to generate charge density
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x <ZrClBr.scf.in >ZrClBr.scf.out

sleep 1

# Step 2: Dense NSCF for DOS
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x <ZrClBr.nscf.in >ZrClBr.nscf.out

sleep 5

# Step 3: Compute total DOS
dos.x <ZrClBr.dos.in >ZrClBr.dos.out

sleep 1

# Step 4: Compute projected DOS (PDOS)
projwfc.x <ZrClBr.pdos.in >ZrClBr.pdos.out
