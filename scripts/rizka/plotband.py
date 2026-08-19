#!/usr/bin/env python3
"""
plot_bands.py — Quantum ESPRESSO Band Structure Plotter
========================================================
Reads Quantum ESPRESSO `bands.x` output and produces a
publication-ready band structure plot.

Usage
-----
    plot_bands.py <prefix> [-o output] [-k labels] [-f Ef] [Emin [Estep [Emax]]]

Positional
----------
    prefix          QE calculation prefix (e.g. "WS2")
                    Looks for  bands.dat.gnu       (band data)
                               <prefix>.band.out   (k-path info)
                               <prefix>.scf.out    (Fermi / VBM energy)
    Emin            Y-axis minimum  (eV, relative to Fermi)   [default: auto]
    Estep           Y-axis tick interval (eV)                  [default: 2]
    Emax            Y-axis maximum  (eV, relative to Fermi)   [default: auto]

Options
-------
    -o, --output    Output filename  (.png / .pdf / .svg)
                    Default: <prefix>_bands.png
    -k, --klabels   Quoted, space-separated high-symmetry labels
                    e.g.  -k "G K M G"   (G is auto-rendered as Gamma)
    -f, --fermi     Override Fermi / VBM energy in eV
    --data          Custom path to bands.dat.gnu
    --dpi           PNG resolution (default: 150)
    --light         Use light (white) theme instead of dark

Examples
--------
    python plot_bands.py WS2
    python plot_bands.py WS2 -o bands.pdf -8 2 2
    python plot_bands.py WS2 -o bands.png -k "G K M G" -8 2 2
    python plot_bands.py WS2 -f -0.6639 -k "G K M G" -o out.svg -10 2 2
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
    band      = "#4fc3f7",
    fermi     = "#ef9a9a",
    fermi_nscf = "#26a69a", # Teal color for NSCF Fermi line
    gap_fill  = "#f9a82522",
    gap_line  = "#f9a825",
    gap_text  = "#f9a825",
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
    band      = "#1565c0",
    fermi     = "#c62828",
    fermi_nscf = "#00796b", # Dark Teal
    gap_fill  = "#fff9c433",
    gap_line  = "#f57f17",
    gap_text  = "#e65100",
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
    """
    Return a "nice" tick step for a given data range.
    Targets approximately n_target intervals, then rounds to the nearest
    clean value: 1, 2, 5, 10, 20, 50, 100 … (powers-of-10 scaled).

    Examples
    --------
    range=78  -> rough=9.75  -> 10
    range=10  -> rough=1.25  -> 2
    range=4   -> rough=0.5   -> 0.5
    range=1.5 -> rough=0.19  -> 0.2
    """
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


# ─── label rendering ────────────────────────────────────────────────────────

_GAMMA_ALIASES = {"g", "gamma", "gam", "г", "γ"}

def _render_label(raw: str) -> str:
    """Convert shorthand k-point names to pretty matplotlib strings."""
    clean = raw.strip()
    if clean.lower() in _GAMMA_ALIASES:
        return r"$\Gamma$"
    # Subscripts: K1 -> K$_1$
    m = re.match(r"^([A-Za-z]+)(\d+)$", clean)
    if m:
        return f"{m.group(1)}$_{{{m.group(2)}}}$"
    return clean


# ─── file parsers ───────────────────────────────────────────────────────────

def _parse_bands_gnu(path: Path):
    """
    Parse bands.dat.gnu -> list of np.ndarray, shape (nk, 2).
    Columns: [k_coordinate, energy_eV]. Bands separated by blank lines.
    """
    bands, buf = [], []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                if buf:
                    bands.append(np.array(buf, dtype=float))
                    buf = []
            else:
                parts = stripped.split()
                if len(parts) == 2:
                    buf.append([float(parts[0]), float(parts[1])])
    if buf:
        bands.append(np.array(buf, dtype=float))
    return bands


def _parse_kpath(band_out: Path):
    """
    Parse high-symmetry points from <prefix>.band.out.
    Returns list of (kx, ky, kz, x_coordinate).
    """
    pattern = re.compile(
        r"high-symmetry point:\s+"
        r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
        r"\s+x coordinate\s+(-?\d+\.\d+)"
    )
    points = []
    if not band_out.exists():
        return points
    with open(band_out) as fh:
        for line in fh:
            m = pattern.search(line)
            if m:
                vals = tuple(float(m.group(i)) for i in range(1, 5))
                points.append(vals)
    return points


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


def _auto_klabels(hs_points):
    """Auto-generate labels: Gamma for (0,0), else indexed."""
    _hex_guesses = ["", "K", "M", "A", "H", "L", "P"]
    labels, counter = [], 1
    for i, (kx, ky, kz, xc) in enumerate(hs_points):
        if abs(kx) < 1e-4 and abs(ky) < 1e-4:
            labels.append(r"$\Gamma$")
        elif len(hs_points) <= len(_hex_guesses) and i < len(_hex_guesses):
            labels.append(_hex_guesses[i])
        else:
            labels.append(f"K{counter}")
            counter += 1
    return labels


def _find_gap(shifted_bands, tol: float = 0.05):
    """
    Auto-detect VBM and CBM directly from shifted band data (VBM = 0 eV).

    Strategy
    --------
    Find the largest band gaps (min of next band - max of current band).
    Among gaps > 0.05 eV, pick the one where the VBM is closest to 0 eV.
    If no gap > 0.05 eV is found, fallback to looking for the highest band
    with max <= tol.

    Returns
    -------
    (vbm_idx, vbm_E, cbm_idx, cbm_E)
    Any of the right-side values can be None if not found.
    """
    band_maxima = [b[:, 1].max() for b in shifted_bands]
    band_minima = [b[:, 1].min() for b in shifted_bands]

    gaps = [band_minima[i+1] - band_maxima[i] for i in range(len(shifted_bands)-1)]
    valid_gaps = [(i, g) for i, g in enumerate(gaps) if g > 0.05]

    if valid_gaps:
        # Sort by how close band_maxima[i] is to 0
        vbm_idx = min(valid_gaps, key=lambda x: abs(band_maxima[x[0]]))[0]
    else:
        # Fallback for metals or no clear gap
        vbm_idx = None
        for i, bmax in enumerate(band_maxima):
            if bmax <= tol:          # band is fully at or below Fermi
                vbm_idx = i
            else:
                break                # bands are ordered; stop at first conduction band

    if vbm_idx is None:
        return None, None, None, None

    vbm_E   = band_maxima[vbm_idx]
    cbm_idx = vbm_idx + 1
    cbm_E   = band_minima[cbm_idx] if cbm_idx < len(shifted_bands) else None
    return vbm_idx, vbm_E, cbm_idx, cbm_E


def _select_centered_bands(shifted_bands, n_show: int, vbm_idx: int):
    """
    Keep n_show bands centred symmetrically around the VBM band.

    Split rule (removes equally from top AND bottom):
        n_val = ceil(n_show / 2)   <- valence bands kept (VBM is the last one)
        n_con = floor(n_show / 2)  <- conduction bands kept (above VBM)

    Example: n_show=12, vbm_idx=12 (0-based), 26 total bands
        n_val=6  -> keep bands 7..12  (indices 7-12)
        n_con=6  -> keep bands 13..18 (indices 13-18)
        total shown = 12
    """
    n_total = len(shifted_bands)
    if n_show <= 0 or n_show >= n_total:
        return shifted_bands, vbm_idx      # no change

    n_val = (n_show + 1) // 2   # valence to keep (incl. VBM)
    n_con = n_show // 2          # conduction to keep

    val_start = max(0, vbm_idx - n_val + 1)
    val_end   = vbm_idx + 1                        # exclusive
    con_end   = min(n_total, vbm_idx + 1 + n_con)

    # If we hit a boundary, compensate on the other side
    if val_start == 0 and (val_end - val_start) < n_val:
        con_end = min(n_total, con_end + (n_val - (val_end - val_start)))
    if con_end == n_total and (con_end - (vbm_idx + 1)) < n_con:
        val_start = max(0, val_start - (n_con - (con_end - (vbm_idx + 1))))

    selected  = shifted_bands[val_start:val_end] + shifted_bands[vbm_idx + 1:con_end]
    new_vbm   = vbm_idx - val_start     # VBM index within the trimmed list
    return selected, new_vbm


# ─── main plotter ───────────────────────────────────────────────────────────

def plot_bands(
    prefix:         str,
    output:         str   = None,
    klabels:        list  = None,
    fermi_in:       float = None,
    data_file:      str   = None,
    Emin:           float = None,
    Estep:          float = None,
    Emax:           float = None,
    dpi:            int   = 150,
    dark:           bool  = False,
    n_bands_show:   int   = None,
    highlight_gap:  bool  = False,
    shift:          bool  = True,
):
    C    = DARK if dark else LIGHT
    cwd  = Path(".")

    # ── locate files ─────────────────────────────────────────────────────
    gnu_path  = Path(data_file) if data_file else cwd / "bands.dat.gnu"
    band_out  = cwd / f"{prefix}.band.out"
    scf_out   = cwd / f"{prefix}.scf.out"
    nscf_out  = cwd / f"{prefix}.nscf.out"
    nscf_dos_out = cwd / f"{prefix}.nscfdos.out"
    out_file  = output or f"{prefix}_bands.png"

    if not gnu_path.exists():
        sys.exit(f"[ERROR] Bands data file not found: {gnu_path}")

    # ── parse bands ──────────────────────────────────────────────────────
    print(f"[INFO] Reading: {gnu_path}")
    bands  = _parse_bands_gnu(gnu_path)
    if not bands:
        sys.exit("[ERROR] No band data parsed.")

    n_bands = len(bands)
    n_kpts  = len(bands[0])
    kx_all  = bands[0][:, 0]
    kmin, kmax = kx_all[0], kx_all[-1]
    print(f"[INFO] Bands: {n_bands},  k-points: {n_kpts}")

    # ── high-symmetry k-path ─────────────────────────────────────────────
    hs_points = _parse_kpath(band_out)
    if not hs_points:
        warnings.warn(f"[WARN] No high-symmetry points in {band_out}")
        hs_x      = [kmin, kmax]
        hs_labels = [r"$\Gamma$", r"$\Gamma$"]
    else:
        hs_x = [p[3] for p in hs_points]
        if klabels:
            diff = len(hs_points) - len(klabels)
            if diff > 0:
                klabels = klabels + ["?"] * diff
            hs_labels = [_render_label(lbl) for lbl in klabels[:len(hs_points)]]
        else:
            hs_labels = _auto_klabels(hs_points)

    print(f"[INFO] k-path: {' -> '.join(hs_labels)}")

    # ── Fermi / VBM energy ───────────────────────────────────────────────
    vbm, cbm, ef, gap_file, is_insulator = _parse_all_energies([scf_out, nscf_dos_out, nscf_out, band_out])
    
    # User override
    if fermi_in is not None:
        vbm = fermi_in
        ef = fermi_in

    # Select reference energy for shifting
    if vbm is not None:
        ref_E = vbm
        ref_type = "VBM"
    elif ef is not None:
        ref_E = ef
        ref_type = "E_F"
    else:
        ref_E = 0.0
        ref_type = "none"

    # Detect gap on VBM-shifted bands (shifted initially by ref_E to find gap near 0)
    vbm_shifted_bands = [b.copy() for b in bands]
    for b in vbm_shifted_bands:
        b[:, 1] -= ref_E
    vbm_idx, vbm_E_shifted, cbm_idx, cbm_E_shifted = _find_gap(vbm_shifted_bands)

    # Calculate true unshifted VBM and CBM from data
    if vbm_idx is not None:
        true_vbm_data = vbm_E_shifted + ref_E
    else:
        true_vbm_data = None

    if cbm_E_shifted is not None:
        true_cbm_data = cbm_E_shifted + ref_E
    else:
        true_cbm_data = None

    # Now define the final shift_val based on true VBM if possible
    if shift:
        if true_vbm_data is not None:
            shift_val = true_vbm_data
            ref_type = "True VBM"
        else:
            shift_val = ref_E
    else:
        shift_val = 0.0

    # Shift bands relative to the final shift_val
    shifted = [b.copy() for b in bands]
    for b in shifted:
        b[:, 1] -= shift_val

    if true_vbm_data is not None:
        vbm_E = true_vbm_data - shift_val
    else:
        vbm_E = None

    if true_cbm_data is not None:
        cbm_E_data = true_cbm_data - shift_val
        eg_data = true_cbm_data - true_vbm_data
        cbm_shifted = cbm_E_data
    else:
        cbm_E_data = None
        eg_data = None
        cbm_shifted = (cbm - shift_val) if (cbm is not None) else None

    if true_vbm_data is not None:
        vbm = true_vbm_data

    print("[INFO] Energy levels parsed from Quantum ESPRESSO output:")
    if vbm is not None:
        print(f"       - Valence Band Maximum (VBM)  : {vbm:+.4f} eV")
    else:
        print("       - Valence Band Maximum (VBM)  : not found in output files")
        
    if ef is not None:
        print(f"       - Fermi Energy (E_F) from NSCF: {ef:+.4f} eV")
    else:
        print("       - Fermi Energy (E_F)          : not found in output files")
        
    if vbm_idx is not None:
        print(f"       - VBM from bands data (0-based idx {vbm_idx}): {vbm_E:+.4f} eV")
    if cbm_E_data is not None:
        print(f"       - CBM from bands data (0-based idx {cbm_idx}): {cbm_E_data:+.4f} eV")
        print(f"       - Band gap from data          : {eg_data:.4f} eV")
    elif gap_file is not None:
        print(f"       - Band gap from outputs       : {gap_file:.4f} eV")

    print(f"[INFO] Shift reference: {ref_type} = {shift_val:+.4f} eV")
    if shift:
        print(f"[INFO] Applied shift: bands shifted by {-shift_val:+.4f} eV")
    else:
        print("[INFO] Applied shift: none (shift disabled)")

    # ── band selection: keep n_bands_show centred at VBM ─────────────────
    if n_bands_show is not None and vbm_idx is not None:
        shifted, vbm_idx = _select_centered_bands(shifted, n_bands_show, vbm_idx)
        n_shown = len(shifted)
        print(f"[INFO] Showing {n_shown} bands (centred at VBM, "
              f"requested {n_bands_show})")
    else:
        n_shown = len(shifted)

    # ── energy window ────────────────────────────────────────────────────
    all_E = np.concatenate([b[:, 1] for b in shifted])
    if Emin is None:
        Emin = np.floor(all_E.min()) - 1.0
    if Emax is None:
        Emax = np.ceil(all_E.max()) + 1.0
    if Estep is None or Estep <= 0:
        Estep = _nice_step(Emax - Emin)   # auto: targets ~8 ticks
    print(f"[INFO] Y range: [{Emin:.1f}, {Emax:.1f}] eV  |  step: {Estep} eV")

    # ── figure ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    fig.patch.set_facecolor(C["fig_bg"])
    ax.set_facecolor(C["axes_bg"])
    for spine in ax.spines.values():
        spine.set_edgecolor(C["border"])
        spine.set_linewidth(1.2)

    # ── band lines ───────────────────────────────────────────────────────
    for i, b in enumerate(shifted):
        ax.plot(b[:, 0], b[:, 1],
                color=C["band"], lw=1.4, alpha=0.90,
                label="Bands" if i == 0 else "_")

    # ── Fermi / VBM horizontal line ──────────────────────────────────────
    if vbm is not None:
        vbm_pos = vbm - shift_val
        vbm_lbl = f"VBM = {vbm_pos:.4f} eV" if not shift else "VBM = 0 eV"
        ax.axhline(vbm_pos, color=C["fermi"], lw=1.6, ls=(0, (8, 4)), label=vbm_lbl, zorder=4)
        
    if ef is not None:
        ef_pos = ef - shift_val
        if not is_insulator and (vbm is None or abs(vbm - ef) > 1e-4):
            # Metal: Ef is physically meaningful, draw it
            ef_lbl = f"$E_F$ = {ef_pos:+.4f} eV" if not shift else f"$E_F$ = {ef_pos:+.4f} eV"
            ax.axhline(ef_pos, color=C["fermi_nscf"], lw=1.6, ls=(0, (3, 3)), label=ef_lbl, zorder=4)
        elif is_insulator:
            print(f"[INFO] Ef line suppressed: insulator detected "
                  f"(Ef = {ef:+.4f} eV from NSCF is a numerical artifact)")

    # ── band gap highlight ────────────────────────────────────────────────
    fermi_y = (vbm - shift_val) if vbm is not None else (ef - shift_val if ef is not None else 0.0)
    gap_E = cbm_shifted
    if highlight_gap and gap_E is None:
        warnings.warn("[WARN] --highlight-gap requested but no CBM found "
                      "(only valence bands in data?). Gap not drawn.")
    if gap_E is not None and Emin < gap_E <= Emax and (highlight_gap or cbm is not None or cbm_E_shifted is not None):
        # Shaded fill between VBM and CBM
        ax.axhspan(fermi_y, gap_E,
                   facecolor=C["gap_fill"], edgecolor="none", zorder=2)
        # CBM dashed line
        ax.axhline(gap_E, color=C["gap_line"], lw=1.5,
                   ls=(0, (4, 3)), zorder=4,
                   label=f"CBM  ($E_g$ = {gap_E - fermi_y:.3f} eV)")
        # Bracket lines: vertical ticks at left edge
        tick_hw = 0.012 * (kmax - kmin)    # half-width of bracket tick
        brak_x  = kmin + 0.055 * (kmax - kmin)
        for yy in (fermi_y, gap_E):
            ax.plot([brak_x - tick_hw, brak_x + tick_hw], [yy, yy],
                    color=C["gap_line"], lw=1.2, zorder=6, solid_capstyle="round")
        # Double-headed arrow
        ax.annotate(
            "", xy=(brak_x, gap_E), xytext=(brak_x, fermi_y),
            arrowprops=dict(arrowstyle="<->", color=C["gap_line"],
                            lw=1.3, shrinkA=0, shrinkB=0),
            zorder=5,
        )
        # Text label
        ax.text(brak_x + 0.022 * (kmax - kmin),
                (fermi_y + gap_E) / 2,
                f"$E_g$ = {gap_E - fermi_y:.3f} eV",
                color=C["gap_text"], va="center", ha="left",
                fontsize=10.5, fontweight="bold", zorder=6)

    # ── high-symmetry lines + labels ────────────────────────────────────
    ax.set_xlim(kmin, kmax)
    ax.set_ylim(Emin, Emax)

    label_y  = Emin - 0.065 * (Emax - Emin)   # below the plot
    for xpos, lbl in zip(hs_x, hs_labels):
        if kmin < xpos < kmax:                   # interior lines only
            ax.axvline(xpos, color=C["border"], lw=1.0, ls="--", zorder=1)
        ax.text(xpos, label_y, lbl,
                ha="center", va="top",
                color=C["text"], fontsize=14, fontweight="bold",
                clip_on=False)

    # ── y-axis ticks ─────────────────────────────────────────────────────
    yticks = np.arange(
        np.ceil(Emin  / Estep) * Estep,
        np.floor(Emax / Estep) * Estep + Estep * 0.01,
        Estep,
    )
    ax.set_yticks(yticks)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.tick_params(axis="y", colors=C["tick"], labelsize=11,
                   which="both", direction="in", right=True)
    ax.set_xticks(hs_x)
    ax.set_xticklabels([""] * len(hs_x))
    ax.tick_params(axis="x", length=0)

    # Grid
    ax.yaxis.grid(True,  color=C["grid"], lw=0.6, zorder=0)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)

    # ── labels / title ───────────────────────────────────────────────────
    ax.set_ylabel("Energy (eV)", color=C["text"], fontsize=13)
    ax.set_title(f"Electronic Band Structure — {prefix}",
                 color=C["title"], fontsize=14, pad=10, fontweight="semibold")

    # ── legend ───────────────────────────────────────────────────────────
    ax.legend(loc="upper right", fontsize=10,
              facecolor=C["legend_bg"], edgecolor=C["border"],
              framealpha=0.85, labelcolor=C["text"])

    # ── info box ─────────────────────────────────────────────────────────
    show_str = f"{n_shown} / {n_bands}" if n_bands_show else str(n_bands)
    info_lines = [
        f"$N_{{\\rm bands}}$ = {show_str}",
        f"$N_k$ = {n_kpts}"
    ]
    if shift:
        ef_str = f"{shift_val:+.4f} eV" if ref_type != "none" else "none"
        info_lines.append(f"Shift = {ef_str}")
    info     = "\n".join(info_lines)
    ax.text(0.99, 0.02, info,
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, color=C["tick"],
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
        prog="plot_bands.py",
        description="Plot Quantum ESPRESSO band structure from bands.dat.gnu",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("prefix",
                   help="QE prefix")
    p.add_argument("-o", "--output", metavar="FILE",
                   help="Output file (.png/.pdf/.svg). Default: <prefix>_bands.png")
    p.add_argument("-k", "--klabels", metavar="LABELS",
                   help='High-symmetry labels in quotes, e.g. "G K M G"')
    p.add_argument("-f", "--fermi", metavar="EV", type=float,
                   help="Override Fermi/VBM energy in eV")
    p.add_argument("--data", metavar="PATH",
                   help="Path to bands.dat.gnu (default: ./bands.dat.gnu)")
    p.add_argument("--dpi", metavar="INT", type=int, default=150,
                   help="PNG DPI (default: 150)")
    p.add_argument("--dark", action="store_true",
                   help="Dark theme")
    p.add_argument("--nbands", metavar="N", type=int, default=None,
                   help=("Show only N bands centred symmetrically around the "
                         "Fermi level. ceil(N/2) valence + floor(N/2) "
                         "conduction bands are kept."))
    p.add_argument("--highlight-gap", action="store_true", dest="highlight_gap",
                   help=("Shade the band gap region and draw a double-arrow "
                         "annotation with Eg. Gap is auto-detected from the "
                         "band data (requires both valence and conduction "
                         "bands to be present)."))
    p.add_argument("--shift", type=str2bool, nargs='?', const=True, default=True,
                   help="Shift energy so Fermi level / VBM is at 0 eV (default: True). Set to False to disable.")
    # Optional positional: Emin  Estep  Emax
    p.add_argument("Emin",  nargs="?", type=float, default=None,
                   help="Y min (eV, relative to Fermi). Default: auto")
    p.add_argument("Estep", nargs="?", type=float, default=None,
                   help="Y tick step (eV). Default: auto (nice round number based on range)")
    p.add_argument("Emax",  nargs="?", type=float, default=None,
                   help="Y max (eV, relative to Fermi). Default: auto")
    return p


def main():
    args    = _build_parser().parse_args()
    klabels = args.klabels.split() if args.klabels else None
    plot_bands(
        prefix         = args.prefix,
        output         = args.output,
        klabels        = klabels,
        fermi_in       = args.fermi,
        data_file      = args.data,
        Emin           = args.Emin,
        Estep          = args.Estep,
        Emax           = args.Emax,
        dpi            = args.dpi,
        dark           = args.dark,
        n_bands_show   = args.nbands,
        highlight_gap  = args.highlight_gap,
        shift          = args.shift,
    )


if __name__ == "__main__":
    main()
