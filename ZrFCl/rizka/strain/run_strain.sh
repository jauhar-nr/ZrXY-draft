#!/bin/bash
#SBATCH --job-name=ZrFCl_STRAIN
#SBATCH --partition=short
#SBATCH --ntasks=32
#SBATCH --output=slurm_strain_%j.log
#SBATCH --error=slurm_strain_%j.err

ulimit -l unlimited
set -e
cd "$SLURM_SUBMIT_DIR"

module load materials/qe/7.2-openmpi

if [ -f "$HOME/.local/miniconda3/etc/profile.d/conda.sh" ]; then 
	source "$HOME/.local/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then 
	source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate qe

# Generate input files jika belum dibuat
python generate_strain_inputs.py

DIRECTIONS=("dir_x" "dir_y")

for d in "${DIRECTIONS[@]}"; do
    for sub in $d/eps_*; do
        if [ -d "$sub" ]; then
            echo "[INFO] Menghitung $sub..."
            pushd "$sub" > /dev/null
            
            mkdir -p tmp
            rm -rf tmp/*
            
            # 1. Run PWscf Relax
            mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < relax.in > relax.out 2>&1
            
            # 2. Run PP (Total Electrostatic Potential untuk Vacuum Level)
            if grep -q "JOB DONE" relax.out; then
                mpirun --use-hwthread-cpus -np $SLURM_NTASKS pp.x < pp.in > pp.out 2>&1
            else
                echo "[WARN] Relax gagal di $sub"
            fi
            
            popd > /dev/null
        fi
    done
done

python analyze_strain.py
