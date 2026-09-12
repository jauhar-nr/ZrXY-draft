#!/bin/bash
#SBATCH --job-name=ZrBr2_BANDS2
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=02:00:00
#SBATCH --output=slurm_bands2_%j.log
#SBATCH --error=slurm_bands2_%j.err

module load gnu12
module load gcc
module load Anaconda
module load materials/qe/7.2-openmpi

echo "============================================================"
echo "  Perhitungan Super-Dense Band Structure: ZrBr2"
echo "  $(date)"
echo "============================================================"

mkdir -p tmp
rm -rf tmp/*

# 1. Eksekusi SCF
echo "[1/4] Menjalankan SCF (32 Cores)..."
mpirun -np 32 pw.x -in scf.in > scf.out
if ! grep -q "JOB DONE" scf.out; then echo "[ERROR] SCF gagal. Berhenti."; exit 1; fi

# 2. Eksekusi NSCFBANDS (Super-Dense, 120 k-points per segment)
echo "[2/4] Menjalankan Super-Dense NSCF Bands (32 Cores, 4 K-point Pools)..."
mpirun -np 32 pw.x -nk 4 -in nscfbands.in > nscfbands.out
if ! grep -q "JOB DONE" nscfbands.out; then echo "[ERROR] NSCFBANDS gagal. Berhenti."; exit 1; fi

# 3. Eksekusi BANDS
echo "[3/4] Menjalankan Ekstraksi BANDS (bands.x)..."
mpirun -np 32 bands.x -in bands.in > bands.out
if ! grep -q "JOB DONE" bands.out; then echo "[ERROR] BANDS gagal. Berhenti."; exit 1; fi

# 4. Ekstraksi Effective Mass & Plotting
echo "[4/4] Menjalankan Analisis Effective Mass & Plotting..."
source ~/.bashrc
conda activate science
python plot_bands.py

echo "============================================================"
echo "  Selesai! Seluruh tahap Super-Dense Band Structure berhasil."
echo "  $(date)"
echo "============================================================"
