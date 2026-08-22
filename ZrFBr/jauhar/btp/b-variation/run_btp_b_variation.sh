#!/bin/bash
#SBATCH --job-name=btp_b_var
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=10:00:00
#SBATCH --output=slurm-%j.out

source ~/.bashrc
conda activate wb

echo "Membuat file interpolasi standar m=7..."
btp2 -vv interpolate -m 7 -e -0.35 -E 0.35 ../tmp/ > /dev/null 2>&1

echo "Menguji variasi b..."
for b in 1000 5000 10000 50000
do
    echo "Running integrate b = $b ..."
    btp2 -vv integrate interpolation.bt2 -b $b 300:1000:100 > log_b_${b}.txt 2>&1
    mv interpolation.trace interpolation_b${b}.trace
done

echo "Membuat plot perbandingan..."
python plot_b_variation.py

echo "Selesai!"
