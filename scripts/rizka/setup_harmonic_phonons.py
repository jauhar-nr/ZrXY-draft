#!/usr/bin/env python3
"""
Setup Harmonic Phonon Calculations (Sub-Fase 2A)
for 2D Janus monolayers: ZrClBr, ZrFBr, and ZrFCl.

Creates 3x3x1 supercell displacements via phonopy-init,
assembles production-ready Quantum ESPRESSO pw.x input files,
and prepares band.conf and Slurm job scripts.
"""

import os
import subprocess
import glob
import re

BASE_DIR = "/mgpfs/home/rabdillah/repo/qe/project/ZrXY/bolz"
PSEUDO_DIR = "/mgpfs/home/rabdillah/repo/qe/pseudo"

MATERIALS = {
    "ZrFBr": {
        "a": 3.344423625,
        "c": 19.812119142,
        "species": [
            ("Zr", 91.224, "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF"),
            ("F", 18.998, "F.pbe-n-kjpaw_psl.1.0.0.UPF"),
            ("Br", 79.904, "Br.pbe-n-kjpaw_psl.1.0.0.UPF")
        ],
        "positions": [
            ("Zr", 0.0000000000, 0.0000000000, 0.5099791267),
            ("F",  0.3333333330, 0.6666666670, 0.5757302267),
            ("Br", 0.3333333330, 0.6666666670, 0.4142906466)
        ]
    },
    "ZrFCl": {
        "a": 3.258871527,
        "c": 18.857449939,
        "species": [
            ("Zr", 91.224, "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF"),
            ("F", 18.998, "F.pbe-n-kjpaw_psl.1.0.0.UPF"),
            ("Cl", 35.453, "Cl.pbe-n-kjpaw_psl.1.0.0.UPF")
        ],
        "positions": [
            ("Zr", 0.0000000000, 0.0000000000, 0.5075110300),
            ("F",  0.3333333330, 0.6666666670, 0.5781975483),
            ("Cl", 0.3333333330, 0.6666666670, 0.4142914217)
        ]
    },
    "ZrClBr": {
        "a": 3.489092374,
        "c": 20.442099183,
        "species": [
            ("Zr", 91.224, "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF"),
            ("Cl", 35.453, "Cl.pbe-n-kjpaw_psl.1.0.0.UPF"),
            ("Br", 79.904, "Br.pbe-n-kjpaw_psl.1.0.0.UPF")
        ],
        "positions": [
            ("Zr", 0.0000000000, 0.0000000000, 0.5027238167),
            ("Cl", 0.3333333330, 0.6666666670, 0.5853014685),
            ("Br", 0.3333333330, 0.6666666670, 0.4119747148)
        ]
    }
}

def create_unitcell_in(mat_name, info, target_dir):
    a = info["a"]
    c = info["c"]
    v1 = [a, 0.0, 0.0]
    v2 = [-a / 2.0, a * (3.0**0.5) / 2.0, 0.0]
    v3 = [0.0, 0.0, c]
    
    lines = [
        "&CONTROL",
        f"  prefix = '{mat_name}',",
        "  outdir = './tmp',",
        f"  pseudo_dir = '{PSEUDO_DIR}',",
        "  calculation = 'scf',",
        "  tprnfor = .TRUE.,",
        "  tstress = .TRUE.,",
        "  disk_io = 'none',",
        "/",
        "&SYSTEM",
        "  ibrav = 0,",
        "  nat   = 3,",
        "  ntyp  = 3,",
        "  ecutwfc = 60,",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-8,",
        "  mixing_beta = 0.7,",
        "/",
        "",
        "ATOMIC_SPECIES"
    ]
    for sp, mass, pseudo in info["species"]:
        lines.append(f"  {sp:<4} {mass:8.3f}   {pseudo}")
    
    lines.append("")
    lines.append("CELL_PARAMETERS (angstrom)")
    lines.append(f"  {v1[0]:14.9f}  {v1[1]:14.9f}  {v1[2]:14.9f}")
    lines.append(f"  {v2[0]:14.9f}  {v2[1]:14.9f}  {v2[2]:14.9f}")
    lines.append(f"  {v3[0]:14.9f}  {v3[1]:14.9f}  {v3[2]:14.9f}")
    lines.append("")
    lines.append("ATOMIC_POSITIONS (crystal)")
    for sp, rx, ry, rz in info["positions"]:
        lines.append(f"  {sp:<4}  {rx:14.10f}  {ry:14.10f}  {rz:14.10f}")
    lines.append("")
    lines.append("K_POINTS {automatic}")
    lines.append("  18 18 1  0 0 0")
    lines.append("")
    
    filepath = os.path.join(target_dir, "unitcell.in")
    with open(filepath, "w") as f:
        f.write("\n".join(lines))
    return filepath

def assemble_supercell_inputs(mat_name, target_dir):
    supercell_files = sorted(glob.glob(os.path.join(target_dir, "supercell-[0-9][0-9][0-9].in")))
    print(f"[{mat_name}] Found {len(supercell_files)} phonopy supercell template files.")
    
    for sc_path in supercell_files:
        basename = os.path.basename(sc_path)
        idx_match = re.search(r'supercell-([0-9]+)\.in', basename)
        if not idx_match:
            continue
        idx_str = idx_match.group(1)
        disp_name = f"disp-{idx_str}"
        
        with open(sc_path, "r") as f:
            sc_content = f.read()
        
        header = [
            "&CONTROL",
            f"  prefix = '{mat_name}_{disp_name}',",
            "  outdir = './tmp',",
            f"  pseudo_dir = '{PSEUDO_DIR}',",
            "  calculation = 'scf',",
            "  tprnfor = .TRUE.,",
            "  tstress = .TRUE.,",
            "  disk_io = 'none',",
            "/",
            "&SYSTEM",
            "  ibrav = 0,",
            "  nat   = 27,",
            "  ntyp  = 3,",
            "  ecutwfc = 60,",
            "/",
            "&ELECTRONS",
            "  conv_thr = 1.0d-8,",
            "  mixing_beta = 0.7,",
            "/",
            ""
        ]
        
        kpoints = [
            "",
            "K_POINTS {automatic}",
            "  6 6 1  0 0 0",
            ""
        ]
        
        full_content = "\n".join(header) + sc_content + "\n".join(kpoints)
        out_in = os.path.join(target_dir, f"{disp_name}.in")
        with open(out_in, "w") as f:
            f.write(full_content)
            
    print(f"[{mat_name}] Generated {len(supercell_files)} production QE inputs (disp-*.in).")

def create_band_conf(mat_name, info, target_dir):
    species_names = " ".join([sp[0] for sp in info["species"]])
    conf = [
        f"ATOM_NAME = {species_names}",
        "DIM = 3 3 1",
        "BAND = 0.0 0.0 0.0  0.5 0.0 0.0  0.3333333333 0.3333333333 0.0  0.0 0.0 0.0",
        "BAND_LABELS = \\Gamma M K \\Gamma",
        "BAND_POINTS = 101",
        "FORCE_CONSTANTS = WRITE"
    ]
    with open(os.path.join(target_dir, "band.conf"), "w") as f:
        f.write("\n".join(conf) + "\n")

def create_local_run_script(mat_name, target_dir):
    content = f"""#!/bin/bash
#SBATCH --job-name=PHONON_{mat_name}
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

echo "=== Memulai Perhitungan Fonon Harmonik: {mat_name} ==="
echo "Host: $(hostname)"
echo "Waktu: $(date)"

mkdir -p tmp
for inp in disp-*.in; do
    out="${{inp%.in}}.out"
    echo ">>> Running pw.x on $inp -> $out ..."
    mpirun --use-hwthread-cpus -np $SLURM_NTASKS pw.x < "$inp" > "$out"
    rm -rf tmp/*
done

echo ">>> Mengumpulkan Forces (phonopy --qe -f)..."
phonopy --qe -f disp-*.out

echo ">>> Menghitung Dispersi Pita Fonon (band.conf)..."
phonopy --qe -c unitcell.in -p band.conf -s

echo "=== Selesai: $(date) ==="
"""
    with open(os.path.join(target_dir, "run_harmonic.job.sh"), "w") as f:
        f.write(content)
    os.chmod(os.path.join(target_dir, "run_harmonic.job.sh"), 0o755)

def main():
    for mat_name, info in MATERIALS.items():
        target_dir = os.path.join(BASE_DIR, mat_name, "phonon", "01_harmonic_phonopy")
        os.makedirs(target_dir, exist_ok=True)
        print(f"\n========================================================")
        print(f"Setting up Harmonic Phonon for {mat_name} in {target_dir}")
        print(f"========================================================")
        
        # 1. Create unitcell.in
        create_unitcell_in(mat_name, info, target_dir)
        
        # 2. Run phonopy-init --qe -d --dim="3 3 1" -c unitcell.in
        cmd = ["phonopy-init", "--qe", "-d", "--dim=3 3 1", "-c", "unitcell.in"]
        print(f"Running: {' '.join(cmd)}")
        res = subprocess.run(cmd, cwd=target_dir, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"[ERROR] phonopy-init failed:\n{res.stderr}")
        else:
            print("[SUCCESS] phonopy-init generated supercells.")
            
        # 3. Assemble full QE disp-*.in
        assemble_supercell_inputs(mat_name, target_dir)
        
        # 4. Create band.conf
        create_band_conf(mat_name, info, target_dir)
        
        # 5. Create local run script
        create_local_run_script(mat_name, target_dir)

    print("\nAll materials successfully set up for Harmonic Phonons!")

if __name__ == "__main__":
    main()
