#!/bin/bash
#SBATCH --job-name=PHONON_ZrFCl
#SBATCH --partition=short
#SBATCH --ntasks=32
#SBATCH --output=slurm_%x_%j.log
#SBATCH --error=slurm_%x_%j.err

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

echo "=== Memulai Perhitungan Fonon Harmonik: ZrFCl ==="
echo "Host: $(hostname)"
echo "Waktu: $(date)"

mkdir -p tmp
for inp in disp-*.in; do
    out="${inp%.in}.out"
    echo ">>> Running pw.x on $inp -> $out ..."
    mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < "$inp" > "$out"
    rm -rf tmp/*
done

echo ">>> Mengumpulkan Forces (phonopy --qe -f)..."
phonopy --qe -f disp-*.out

echo ">>> Menghitung Dispersi Pita Fonon (band.conf)..."
phonopy --qe -c unitcell.in -p band.conf -s

echo "=== Selesai: $(date) ==="
