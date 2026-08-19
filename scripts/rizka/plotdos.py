#!/usr/bin/env python3
"""
plotdos.py — Quantum ESPRESSO Density of States (DOS) Plotter
==============================================================
Reads Quantum ESPRESSO `dos.x` output and produces a density of states plot.

Usage
-----
    plotdos.py <prefix> [-o output] [-f Ef] [--data path] [Emin [Estep [Emax]]]

Positional
----------
    prefix          QE calculation prefix (e.g. "WS2_1")
                    Looks for  dos.dat             (DOS data)
                               <prefix>.nscfdos.out (Fermi / VBM energy)
                               <prefix>.scf.out    (Fermi / VBM energy)
    Emin            Energy axis minimum (eV, relative to Fermi)   [default: auto]
    Estep           Energy axis tick interval (eV)                 [default: 2]
    Emax            Energy axis maximum (eV, relative to Fermi)   [default: auto]

Options
-------
    -o, --output    Output filename (.png / .pdf / .svg)
                    Default: <prefix>_dos.png
    -f, --fermi     Override Fermi / VBM energy in eV
    --data          Custom path to DOS data file (default: ./dos.dat)
    --dpi           PNG resolution (default: 150)
    --light         Use light (white) theme instead of dark
    --energy-y      Plot with Energy on the Y-axis (aligns with band structure plots)
    --no-fill       Do not shade the area under the DOS curve
"""

import sys
import re
import argparse
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ─── colour palettes ────────────────────────────────────────────────────────

DARK = dict(
    fig_bg    = "#0f1117",
    axes_bg   = "#161b27",
    border    = "#334155",
    dos_line  = "#4fc3f7",  # Matching band color in plot_bands.py
    dos_fill  = "#4fc3f722", # Shading under DOS curve
    fermi     = "#ef9a9a",  # Matching Fermi line in plot_bands.py
    fermi_nscf = "#26a69a", # Teal color for NSCF Fermi line
    text      = "#e2e8f0",
    tick      = "#94a3b8",
    grid      = "#1e2a3a",
    legend_bg = "#1e2a3a",
    title     = "#e2e8f0",
)

LIGHT = dict(
    fig_bg    = "#f8fafc",
    axes_bg   = "#ffffff",
    border    = "#64748b",
    dos_line  = "#1565c0",
    dos_fill  = "#1565c022",
    fermi     = "#c62828",
    fermi_nscf = "#00796b", # Dark Teal
    text      = "#1e293b",
    tick      = "#475569",
    grid      = "#e2e8f0",
    legend_bg = "#f1f5f9",
    title     = "#0f172a",
)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# ─── tick helpers ──────────────────────────────────────────────────────────

def _nice_step(data_range: float, n_target: int = 8) -> float:
    """Return a 'nice' tick step for a given data range."""
    if data_range <= 0:
        return 1.0
    rough = data_range / n_target
    mag   = 10.0 ** np.floor(np.log10(rough))
    frac  = rough / mag
    if   frac <= 1.0: nice = 1.0
    elif frac <= 2.0: nice = 2.0
    elif frac <= 5.0: nice = 5.0
    else:             nice = 10.0
    return nice * mag


# ─── file parsers ───────────────────────────────────────────────────────────

def _parse_dos_file(path: Path):
    """
    Parse DOS data file.
    Expected columns: [Energy_eV, DOS_states_per_eV, Integrated_DOS]
    """
    data = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    data.append([float(x) for x in parts[:3]])
                except ValueError:
                    continue
    return np.array(data)


def _parse_all_energies(file_list: list):
    """
    Search all QE output files to find VBM, CBM, and Fermi energy (Ef).
    Returns (vbm, cbm, ef, gap_from_file, is_insulator).
    is_insulator is True when QE printed 'highest occupied' (insulator/semiconductor).
    """
    vbm, cbm, ef, gap_file, is_insulator = None, None, None, None, False
    for fpath in file_list:
        if not fpath or not fpath.exists():
            continue
        text = fpath.read_text(errors="replace")
        
        # 1. Look for VBM and CBM printed together:
        m = re.search(
            r"highest occupied,\s*lowest unoccupied level.*?([-\d.]+)\s+([-\d.]+)",
            text, re.I
        )
        if m:
            vbm_val, cbm_val = float(m.group(1)), float(m.group(2))
            if vbm is None: vbm = vbm_val
            if cbm is None: cbm = cbm_val
            is_insulator = True
            
        # 2. Look for VBM alone:
        m = re.search(r"highest occupied level\s*\(ev\)\s*:\s*([-\d.]+)", text, re.I)
        if m:
            vbm_val = float(m.group(1))
            if vbm is None: vbm = vbm_val
            is_insulator = True

        # 3. Look for Fermi energy:
        m = re.search(r"the Fermi energy is\s+([-\d.]+)\s+ev", text, re.I)
        if m:
            ef_val = float(m.group(1))
            if ef is None: ef = ef_val
        m = re.search(r"Fermi energy\s*=\s*([-\d.]+)\s+ev", text, re.I)
        if m:
            ef_val = float(m.group(1))
            if ef is None: ef = ef_val
            
    if vbm is not None and cbm is not None:
        gap_file = cbm - vbm
        
    return vbm, cbm, ef, gap_file, is_insulator


# ─── main plotter ───────────────────────────────────────────────────────────

def plot_dos(
    prefix:        str,
    output:        str   = None,
    fermi_in:      float = None,
    data_file:     str   = None,
    Emin:          float = None,
    Estep:         float = None,
    Emax:          float = None,
    dpi:           int   = 150,
    dark:          bool  = False,
    energy_y:      bool  = False,
    no_fill:       bool  = False,
    shift:         bool  = True,
):
    C    = DARK if dark else LIGHT
    cwd  = Path(".")

    # ── locate files ─────────────────────────────────────────────────────
    data_path = Path(data_file) if data_file else cwd / "dos.dat"
    band_out  = cwd / f"{prefix}.band.out"
    scf_out   = cwd / f"{prefix}.scf.out"
    nscf_out  = cwd / f"{prefix}.nscfdos.out"
    # Fallback to standard nscf.out if nscfdos.out is not found
    if not nscf_out.exists():
        nscf_out = cwd / f"{prefix}.nscf.out"
        
    out_file  = output or f"{prefix}_dos.png"

    if not data_path.exists():
        sys.exit(f"[ERROR] DOS data file not found: {data_path}")

    # ── parse DOS data ───────────────────────────────────────────────────
    print(f"[INFO] Reading: {data_path}")
    raw_data = _parse_dos_file(data_path)
    if len(raw_data) == 0:
        sys.exit("[ERROR] No DOS data parsed.")

    # ── Fermi / VBM energy ───────────────────────────────────────────────
    nscf_dos_out = cwd / f"{prefix}.nscfdos.out"
    vbm, cbm, ef, gap_file, is_insulator = _parse_all_energies([scf_out, nscf_dos_out, nscf_out, band_out])
    
    # User override
    if fermi_in is not None:
        vbm = fermi_in
        ef = fermi_in

    # Select reference energy for shifting (prefer VBM)
    if vbm is not None:
        ref_E = vbm
        ref_type = "VBM"
    elif ef is not None:
        ref_E = ef
        ref_type = "E_F"
    else:
        ref_E = 0.0
        ref_type = "none"

    # Define actual shift value
    shift_val = ref_E if shift else 0.0

    print("[INFO] Energy levels parsed from Quantum ESPRESSO output:")
    if vbm is not None:
        print(f"       - Valence Band Maximum (VBM)  : {vbm:+.4f} eV")
    else:
        print("       - Valence Band Maximum (VBM)  : not found in output files")
        
    if ef is not None:
        print(f"       - Fermi Energy (E_F) from NSCF: {ef:+.4f} eV")
    else:
        print("       - Fermi Energy (E_F)          : not found in output files")
        
    if cbm is not None:
        print(f"       - Conduction Band Minimum (CBM): {cbm:+.4f} eV")
        
    if gap_file is not None:
        print(f"       - Band Gap (from outputs)     : {gap_file:.4f} eV")
    if is_insulator:
        print(f"       - System type                 : insulator/semiconductor")
        if ef is not None:
            print(f"       - NOTE: Ef = {ef:+.4f} eV is a tetrahedra artifact (suppressed from plot)")
        
    print(f"[INFO] Shift reference: {ref_type} = {ref_E:+.4f} eV")
    if shift:
        print(f"[INFO] Applied shift: energies shifted by {-shift_val:+.4f} eV")
    else:
        print("[INFO] Applied shift: none (shift disabled)")

    # Shift energies relative to Fermi level
    energies_shifted = raw_data[:, 0] - shift_val
    dos_vals = raw_data[:, 1]

    # ── energy window ────────────────────────────────────────────────────
    if Emin is None:
        Emin = np.floor(energies_shifted.min()) - 1.0
    if Emax is None:
        Emax = np.ceil(energies_shifted.max()) + 1.0
    if Estep is None or Estep <= 0:
        Estep = _nice_step(Emax - Emin)
    print(f"[INFO] Y range: [{Emin:.1f}, {Emax:.1f}] eV  |  step: {Estep} eV")

    # ── figure ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    fig.patch.set_facecolor(C["fig_bg"])
    ax.set_facecolor(C["axes_bg"])
    for spine in ax.spines.values():
        spine.set_edgecolor(C["border"])
        spine.set_linewidth(1.2)

    # ── plotting curve ───────────────────────────────────────────────────
    if energy_y:
        # Energy on Y-axis, DOS on X-axis (aligns vertically with bands)
        ax.plot(dos_vals, energies_shifted, color=C["dos_line"], lw=1.6, label="DOS")
        if not no_fill:
            ax.fill_betweenx(energies_shifted, 0, dos_vals, facecolor=C["dos_fill"], edgecolor="none")
        
        # Fermi / VBM reference lines
        if vbm is not None:
            vbm_pos = vbm - shift_val
            vbm_lbl = f"VBM = {vbm_pos:.4f} eV" if not shift else "VBM = 0 eV"
            ax.axhline(vbm_pos, color=C["fermi"], lw=1.6, ls=(0, (8, 4)), label=vbm_lbl, zorder=4)
        if ef is not None:
            ef_pos = ef - shift_val
            if not is_insulator and (vbm is None or abs(vbm - ef) > 1e-4):
                ef_lbl = f"$E_F$ = {ef_pos:+.4f} eV" if not shift else f"$E_F$ = {ef_pos:+.4f} eV"
                ax.axhline(ef_pos, color=C["fermi_nscf"], lw=1.6, ls=(0, (3, 3)), label=ef_lbl, zorder=4)

        # Set axes labels and limits
        ax.set_ylabel("Energy (eV)", color=C["text"], fontsize=13)
        ax.set_xlabel("Density of States (states/eV)", color=C["text"], fontsize=13)
        ax.set_ylim(Emin, Emax)
        ax.set_xlim(0, dos_vals[(energies_shifted >= Emin) & (energies_shifted <= Emax)].max() * 1.05)
        
        # Ticks formatting
        yticks = np.arange(
            np.ceil(Emin  / Estep) * Estep,
            np.floor(Emax / Estep) * Estep + Estep * 0.01,
            Estep,
        )
        ax.set_yticks(yticks)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.tick_params(axis="both", colors=C["tick"], labelsize=11, which="both", direction="in", right=True, top=True)
        ax.yaxis.grid(True, color=C["grid"], lw=0.6, zorder=0)
        
    else:
        # Energy on X-axis, DOS on Y-axis (standard plot)
        ax.plot(energies_shifted, dos_vals, color=C["dos_line"], lw=1.6, label="DOS")
        if not no_fill:
            ax.fill_between(energies_shifted, 0, dos_vals, facecolor=C["dos_fill"], edgecolor="none")
        
        # Fermi / VBM reference lines
        if vbm is not None:
            vbm_pos = vbm - shift_val
            vbm_lbl = f"VBM = {vbm_pos:.4f} eV" if not shift else "VBM = 0 eV"
            ax.axvline(vbm_pos, color=C["fermi"], lw=1.6, ls=(0, (8, 4)), label=vbm_lbl, zorder=4)
        if ef is not None:
            ef_pos = ef - shift_val
            if not is_insulator and (vbm is None or abs(vbm - ef) > 1e-4):
                ef_lbl = f"$E_F$ = {ef_pos:+.4f} eV" if not shift else f"$E_F$ = {ef_pos:+.4f} eV"
                ax.axvline(ef_pos, color=C["fermi_nscf"], lw=1.6, ls=(0, (3, 3)), label=ef_lbl, zorder=4)

        # Set axes labels and limits
        ax.set_xlabel("Energy (eV)", color=C["text"], fontsize=13)
        ax.set_ylabel("Density of States (states/eV)", color=C["text"], fontsize=13)
        ax.set_xlim(Emin, Emax)
        ax.set_ylim(0, dos_vals[(energies_shifted >= Emin) & (energies_shifted <= Emax)].max() * 1.05)

        # Ticks formatting
        xticks = np.arange(
            np.ceil(Emin  / Estep) * Estep,
            np.floor(Emax / Estep) * Estep + Estep * 0.01,
            Estep,
        )
        ax.set_xticks(xticks)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.tick_params(axis="both", colors=C["tick"], labelsize=11, which="both", direction="in", right=True, top=True)
        ax.xaxis.grid(True, color=C["grid"], lw=0.6, zorder=0)

    ax.set_axisbelow(True)

    # ── labels / title ───────────────────────────────────────────────────
    ax.set_title(f"Density of States — {prefix}", color=C["title"], fontsize=14, pad=10, fontweight="semibold")

    # ── legend ───────────────────────────────────────────────────────────
    ax.legend(loc="upper right", fontsize=10,
              facecolor=C["legend_bg"], edgecolor=C["border"],
              framealpha=0.85, labelcolor=C["text"])

    # ── info box ─────────────────────────────────────────────────────────
    if shift:
        ef_str = f"{ref_E:+.4f} eV" if ref_type != "none" else "none"
        info = f"Shift = {ef_str}"
        ax.text(0.98, 0.02, info,
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=9.5, color=C["tick"],
                bbox=dict(boxstyle="round,pad=0.4", facecolor=C["legend_bg"],
                          edgecolor=C["border"], alpha=0.8))

    # ── save ─────────────────────────────────────────────────────────────
    plt.tight_layout(pad=1.2)
    suffix  = Path(out_file).suffix.lower()
    save_kw = dict(dpi=dpi, bbox_inches="tight", facecolor=C["fig_bg"])
    if suffix in (".pdf", ".svg", ".eps"):
        save_kw.pop("dpi", None)

    fig.savefig(out_file, **save_kw)
    plt.close(fig)
    print(f"[OK]   Saved -> {out_file}")


# ─── CLI ────────────────────────────────────────────────────────────────────

def _build_parser():
    p = argparse.ArgumentParser(
        prog="plotdos.py",
        description="Plot Quantum ESPRESSO density of states (DOS) from dos.dat",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("prefix",
                   help="QE prefix")
    p.add_argument("-o", "--output", metavar="FILE",
                   help="Output file (.png/.pdf/.svg). Default: <prefix>_dos.png")
    p.add_argument("-f", "--fermi", metavar="EV", type=float,
                   help="Override Fermi/VBM energy in eV")
    p.add_argument("--data", metavar="PATH",
                   help="Path to DOS data file (default: ./dos.dat)")
    p.add_argument("--dpi", metavar="INT", type=int, default=150,
                   help="PNG DPI (default: 150)")
    p.add_argument("--dark", action="store_true",
                   help="Dark theme")
    p.add_argument("--energy-y", action="store_true",
                   help="Plot energy on the Y-axis (ideal for side-by-side band comparison)")
    p.add_argument("--no-fill", action="store_true",
                   help="Do not fill/shade the area under the DOS curve")
    p.add_argument("--shift", type=str2bool, nargs='?', const=True, default=True,
                   help="Shift energy so Fermi level / VBM is at 0 eV (default: True). Set to False to disable.")
    # Optional positional: Emin  Estep  Emax
    p.add_argument("Emin",  nargs="?", type=float, default=None,
                   help="Energy min (relative to Fermi). Default: auto")
    p.add_argument("Estep", nargs="?", type=float, default=None,
                   help="Energy tick step. Default: auto")
    p.add_argument("Emax",  nargs="?", type=float, default=None,
                   help="Energy max (relative to Fermi). Default: auto")
    return p


def main():
    args = _build_parser().parse_args()
    plot_dos(
        prefix    = args.prefix,
        output    = args.output,
        fermi_in  = args.fermi,
        data_file = args.data,
        Emin      = args.Emin,
        Estep     = args.Estep,
        Emax      = args.Emax,
        dpi       = args.dpi,
        dark      = args.dark,
        energy_y  = args.energy_y,
        no_fill   = args.no_fill,
        shift     = args.shift,
    )


if __name__ == "__main__":
    main()
