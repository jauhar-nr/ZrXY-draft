#!/bin/bash
#SBATCH --job-name=btp_m_var
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=10:00:00
#SBATCH --output=slurm-%j.out

source ~/.bashrc
conda activate wb

echo "Menguji variasi m..."
rm -f summary_m.txt

for m in {2..32}
do
    echo "Running interpolate m = $m ..."
    btp2 -vv interpolate -m $m -e -0.35 -E 0.35 ../tmp/ > log_m_${m}.txt 2>&1
    
    kpts=$(grep "irreducible k points have been generated" log_m_${m}.txt | awk '{print $3}')
    echo "Jumlah k-points untuk m=$m : $kpts" >> summary_m.txt
    
    rm -f *.bt2
done
echo "Selesai!"
