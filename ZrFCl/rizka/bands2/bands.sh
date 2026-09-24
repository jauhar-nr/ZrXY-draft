#!/bin/bash
#SBATCH --job-name=ZrFCl_BANDS2
#SBATCH --partition=short
#SBATCH --ntasks=32
#SBATCH --output=slurm_bands2_%j.log
#SBATCH --error=slurm_bands2_%j.err

ulimit -l unlimited
set -e
cd "$SLURM_SUBMIT_DIR"

module load materials/qe/7.2-openmpi

mkdir -p tmp
rm -rf tmp/*

# 1. Eksekusi SCF
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < scf.in > scf.out

# 2. Eksekusi NSCFBANDS (Super-Dense, 120 k-points per segment)
mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x -nk 4 < nscfbands.in > nscfbands.out

# 3. Eksekusi BANDS
mpirun --use-hwthread-cpus -np $SLURM_NTASKS bands.x < bands.in > bands.out

# 4. Ekstraksi Effective Mass & Plotting
if [ -f "$HOME/.local/miniconda3/etc/profile.d/conda.sh" ]; then 
	source "$HOME/.local/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then 
	source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate qe
python plot_bands.py
