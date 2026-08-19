#!/bin/bash

#SBATCH --job-name="ZrFCL dos+pdos"
#SBATCH --partition=short
#SBATCH --ntasks=32
#SBATCH --output=slurm_ZrFCl_new_dos%j.log
#SBATCH --error=slurm_ZrFCl_new_dos%j.err

# NOTE: Runs SCF then dense NSCF, then computes total DOS and projected DOS (PDOS).

ulimit -l unlimited

module load materials/qe/7.2-openmpi

# Go to the dos directory
cd $SLURM_SUBMIT_DIR

# Step 1: SCF to generate charge density
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x <ZrFCl.scf.in >ZrFCl.scf.out

sleep 5

# Step 2: Dense NSCF for DOS
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x <ZrFCl.nscf.in >ZrFCl.nscf.out

sleep 5

# Step 3: Compute total DOS
dos.x <ZrFCl.dos.in >ZrFCl.dos.out

sleep 5

# Step 4: Compute projected DOS (PDOS)
projwfc.x <ZrFCl.pdos.in >ZrFCl.pdos.out
