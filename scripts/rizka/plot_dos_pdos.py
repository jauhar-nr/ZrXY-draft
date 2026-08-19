#!/bin/python
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scienceplots

plt.style.use(['science', 'no-latex']) # Use science style without requiring LaTeX installation

#
#
import os

# def get_efermi(_file):
#     if not os.path.exists(_file):
#         raise FileNotFoundError(f"File not found: {scf_file}")
#
#     feermi = None
#     with open(scf_file, 'r') as file:
#         for line in file:
#             if "the Fermi energy is" in line:
#                 fermi_energy = float(line.split()[4])
#
#     if fermi_energy is None:
#         print("Warning: Fermi energy not found.")
#
#     return fermi_energy
#
# def get_occupy(nscf_file):
#     if not os.path.exists(nscf_file):
#         raise FileNotFoundError(f"File not found: {nscf_file}")
#
#     vbm, cbm = None, None
#     with open(nscf_file, 'r') as file:
#         for line in file:
#             if "highest occupied, lowest unoccupied level (ev):" in line:
#                 values = line.split(":")[1].split()
#                 vbm = float(values[0])
#                 cbm = float(values[1])
#
#             elif "highest occupied level (ev):" in line:
#                 values = line.split(":")[1].split()
#                 vbm = float(values[0])
#                 print("Only VBM (highest occupied) found.")
#
#     if vbm is None and cbm is None:
#         print("Warning: VBM/CBM not found.")
#
#     return vbm, cbm
#
# ── colour palette ─────────────────────────────────────
C = {
    "total": "#b0b0b0",
    "Zr_s":  "#90caf9",   # light blue
    "Zr_p":  "#42a5f5",   # blue
    "Zr_d":  "#0d47a1",   # deep blue
    "X_s":   "#ffe082",   # amber
    "X_p":   "#e65100",   # deep orange
    "Y_s":   "#a5d6a7",   # light green
    "Y_p":   "#1b5e20",   # dark green
}

MATERIALS = [
    {
        "label":  "ZrClBr",
        "prefix": "ZrClBr",
        "X":      "Cl",
        "Y":      "Br",
        "efermi": 5.7994,   # dos.dat EFermi (= VBM, NOT mid-gap)
        "vbm":    4.7913,   # scf.out highest occupied level (eV)
        "cbm":    5.8064,   # scf.out lowest unoccupied (< VBM — coarse 10×10 SCF artifact)
        "path":   Path("ZrClBr/dos/new_dos"),
    },
    {
        "label":  "ZrFBr",
        "prefix": "ZrFBr",
        "X":      "F",
        "Y":      "Br",
        "efermi": 5.5630,   # dos.dat EFermi (= VBM, NOT mid-gap)
        "vbm":    4.2921,   # scf.out highest occupied level (eV)
        "cbm":    5.5808,   # scf.out lowest unoccupied (< VBM — coarse 10×10 SCF artifact)
        "path":   Path("ZrFBr/dos/new_dos"),
    },
    {
        "label":  "ZrFCl",
        "prefix": "ZrFCl",
        "X":      "F",
        "Y":      "Cl",
        "efermi": 5.2926,   # dos.dat EFermi (= VBM, NOT mid-gap)
        "vbm":    3.9298,   # scf.out highest occupied level (eV)
        "cbm":    5.2920,   # scf.out lowest unoccupied (< VBM — coarse 10×10 SCF artifact)
        "path":   Path("ZrFCl/dos/new_dos"),
    },
]

EMIN, EMAX = -7.0, 4.0   # relative to E_F

# plt.rcParams.update({
#     "font.family":      "DejaVu Sans",
#     "font.size":        11,
#     "axes.linewidth":   1.2,
#     "axes.spines.top":  False,
#     "axes.spines.right":False,
#     "xtick.direction":  "in",
#     "ytick.direction":  "in",
#     "xtick.major.size": 4,
#     "ytick.major.size": 4,
# })

# ── func ───────────────────────────────────────────────────────────────────
def load_ldos(files):
    """Sum ldos column (col 1) across given file list."""
    data = None
    for f in files:
        d = np.loadtxt(f, comments="#")
        data = d[:, 1] if data is None else data + d[:, 1]
    return data

# ── one figure per material ───────────────────────────────────────────────────
for mat in MATERIALS:
    p       = mat["path"]
    ef      = mat["efermi"]
    vbm_rel = mat["vbm"] - ef   # VBM relative to E_F
    cbm_rel = mat["cbm"] - ef   # CBM relative to E_F
    X, Y    = mat["X"], mat["Y"]
    label   = mat["label"]

    # energy grid from pdos_tot (1 row more than dos.dat — use as master)
    pdos_tot_data = np.loadtxt(p / "pdos.pdos_tot", comments="#")
    E_raw = pdos_tot_data[:, 0]
    E     = E_raw - ef
    mask  = (E >= EMIN) & (E <= EMAX)
    E     = E[mask]
    # total PDOS sum as reported by projwfc.x (col 1 of pdos_tot)
    pdos_tot_sum = pdos_tot_data[:, 1][mask]

    def ldos(fnames):
        return load_ldos([p / f for f in fnames])[mask]

    # total DOS — interpolate dos.dat onto pdos energy grid
    dos_data  = np.loadtxt(p / "dos.dat", comments="#")
    dos_total = np.interp(E_raw[mask], dos_data[:, 0], dos_data[:, 1])

    # Zr
    zr_s = ldos([f"pdos.pdos_atm#1(Zr)_wfc#1(s)",
                  f"pdos.pdos_atm#1(Zr)_wfc#2(s)"])
    zr_p = ldos([f"pdos.pdos_atm#1(Zr)_wfc#3(p)",
                  f"pdos.pdos_atm#1(Zr)_wfc#4(p)"])
    zr_d = ldos([f"pdos.pdos_atm#1(Zr)_wfc#5(d)"])

    # halogen X (top, atom #2)
    x_s = ldos([f"pdos.pdos_atm#2({X})_wfc#1(s)"])
    x_p = ldos([f"pdos.pdos_atm#2({X})_wfc#2(p)"])

    # halogen Y (bottom, atom #3)
    y_s = ldos([f"pdos.pdos_atm#3({Y})_wfc#1(s)"])
    y_p = ldos([f"pdos.pdos_atm#3({Y})_wfc#2(p)"])

    # auto y limit: use 97th-percentile to avoid isolated VHS spikes
    ymax = np.ceil(np.percentile(dos_total, 97) / 2) * 2 + 1

    # ── figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # total DOS shaded — from dos.x (plane-wave basis, all states)
    ax.fill_between(E, dos_total, 0, alpha=0.25, color=C["total"], linewidth=0)
    ax.plot(E, dos_total, color="#888888", lw=0.9, label="Total (dos.x)")

    # projwfc.x projection sum — should match Total in valence, less so in conduction
    # ax.plot(E, pdos_tot_sum, color="#111111", lw=0.8, ls="--",
            # label="Proj. sum (projwfc)", alpha=0.7)

    # PDOS lines — (energy, dos)
    ax.plot(E, zr_d, color=C["Zr_d"], lw=1.8, label="Zr $d$")
    ax.plot(E, zr_p, color=C["Zr_p"], lw=1.4, label="Zr $p$")
    ax.plot(E, zr_s, color=C["Zr_s"], lw=1.1, label="Zr $s$", ls="--")
    ax.plot(E, x_p,  color=C["X_p"],  lw=1.8, label=f"{X} $p$")
    ax.plot(E, x_s,  color=C["X_s"],  lw=1.1, label=f"{X} $s$", ls="--")
    ax.plot(E, y_p,  color=C["Y_p"],  lw=1.8, label=f"{Y} $p$")
    ax.plot(E, y_s,  color=C["Y_s"],  lw=1.1, label=f"{Y} $s$", ls="--")

    # Fermi level — vertical line at E=0
    ax.axvline(0, color="#d32f2f", lw=1.0, ls="--", alpha=0.8, label=f"$E_F$ = {ef:.3f} eV")

    # VBM marker
    ax.axvline(vbm_rel, color="#6a1b9a", lw=1.0, ls=":", alpha=0.9,
               label=f"VBM = {mat['vbm']:.4f} eV ({vbm_rel:+.4f})")

    # CBM marker
    ax.axvline(cbm_rel, color="#2e7d32", lw=1.0, ls=":", alpha=0.9,
               label=f"CBM = {mat['cbm']:.4f} eV ({cbm_rel:+.4f})")

    # axes — E on X, DOS on Y
    ax.set_xlim(EMIN, EMAX)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("$E - E_F$ (eV)", fontsize=12)
    ax.set_ylabel("DOS (states/eV)", fontsize=12)
    ax.set_title(label, fontsize=14, fontweight="bold", pad=10)

    # x ticks every 2 eV along energy axis
    ax.set_xticks(range(int(EMIN), int(EMAX) + 1, 1))

    # legend
    ax.legend(loc="upper right", fontsize=9, framealpha=1,
              facecolor="white", edgecolor="#cccccc", bbox_to_anchor=(1, 1))

    # grid (subtle vertical lines along energy axis)
    ax.xaxis.grid(True, color="#e0e0e0", lw=0.6, ls="-")
    ax.set_axisbelow(True)

    ax.tick_params(labelsize=10)

    fig.tight_layout()
    out = Path(f"dos_pdos_{label}w.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved → {out.resolve()}")
    plt.close(fig)

print("All done.")
