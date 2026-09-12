#!/bin/bash
#SBATCH --job-name=ZrBr2_STRAIN
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=04:00:00
#SBATCH --output=slurm_strain_%j.log
#SBATCH --error=slurm_strain_%j.err

module load gnu12
module load gcc
module load Anaconda
module load materials/qe/7.2-openmpi

echo "============================================================"
echo "  Mulai Perhitungan Regangan Kisi (DPT): ZrBr2"
echo "  $(date)"
echo "============================================================"

DIRECTIONS=("dir_x" "dir_y")

for d in "${DIRECTIONS[@]}"; do
    echo "------------------------------------------------------------"
    echo "  Memproses Arah: $d"
    echo "------------------------------------------------------------"
    
    for sub in $d/eps_*; do
        if [ -d "$sub" ]; then
            echo "[INFO] Menghitung $sub..."
            pushd "$sub" > /dev/null
            
            mkdir -p tmp
            rm -rf tmp/*
            
            # 1. Run PWscf Relax
            mpirun -np 32 pw.x -in relax.in > relax.out 2>&1
            
            # 2. Run PP (Total Electrostatic Potential untuk Vacuum Level)
            if grep -q "JOB DONE" relax.out; then
                mpirun -np 32 pp.x -in pp.in > pp.out 2>&1
            else
                echo "[WARN] Relax gagal di $sub"
            fi
            
            popd > /dev/null
        fi
    done
done

echo "------------------------------------------------------------"
echo "  Menganalisis Hasil Regangan (C, E_1, tau, mu)..."
echo "------------------------------------------------------------"
source ~/.bashrc
conda activate science
python analyze_strain.py

echo "============================================================"
echo "  Selesai Seluruh Perhitungan Regangan ZrBr2!"
echo "  $(date)"
echo "============================================================"
