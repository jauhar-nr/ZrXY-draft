#!/bin/bash
#SBATCH --job-name=btp_m_plot
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=10:00:00
#SBATCH --output=slurm_m_plot-%j.out

source ~/.bashrc
conda activate wb

echo "Mulai menguji integrate dan plot untuk m=2, 16, dan 32..."

for m in 2 16 32
do
    echo "Running interpolate and integrate m = $m ..."
    btp2 -vv interpolate -m $m -e -0.35 -E 0.35 ../tmp/ > /dev/null 2>&1
    
    # Jalankan integrate menggunakan konfigurasi run_hpc_btp.sh (-b 10000 300:1000:300)
    btp2 -vv integrate interpolation.bt2 -b 10000 300:1000:300 > /dev/null 2>&1
    
    mv interpolation.trace interpolation_m${m}.trace
    rm -f *.bt2
done

echo "Membuat plot 3-panel..."
python plot_m_variation.py

echo "Selesai!"
