#!/bin/bash
#SBATCH --job-name=btp_ZrFBr
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# 1. Mengikuti workflow HPC Anda
module load gnu12
module load gcc
module load materials/qe/7.2-openmpi

echo "Mulai Kalkulasi SCF..."
mpirun -np 32 pw.x -in scf.in > scf.out

echo "Mulai Kalkulasi NSCF..."
mpirun -np 32 pw.x -in nscf.in > nscf.out

# Pindah sementara ke ekosistem python untuk BoltzTraP2 (Opsional, tergantung setting HPC)
# btp2 biasanya ada di env python Anda, jadi kita harus mengaktifkannya
echo "Mulai Analisis BoltzTraP2..."
source ~/.bashrc
conda activate wb

# Hapus folder BTP lama jika ada agar bersih
rm -rf ./tmp/*.bt2

# Interpolasi BoltzTraP2 menggunakan multiplier 7
btp2 -vv interpolate -m 7 -e -0.35 -E 0.35 ./tmp/

# Integrasi dari 300K ke 1000K (agar 900K masuk)
btp2 -vv integrate interpolation.bt2 -b 10000 300:1000:300

# Plot
python plot_btp.py

echo "Selesai Keseluruhan Workflow HPC!"
