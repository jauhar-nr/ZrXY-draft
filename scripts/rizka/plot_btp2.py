#!/usr/bin/env python3
"""
plot_btp2.py  —  Plot BoltzTraP2 thermoelectric properties from interpolation.condtens

Usage:
    python plot_btp2.py PREFIX E_FERMI [options]

Positional arguments:
    PREFIX        Material prefix (e.g. ZrClBr). Used for output filenames and title.
    E_FERMI       Fermi / midgap energy in eV  [= (VBM + CBM) / 2  from *.scf.out]
                  Ignored if --trace is given (refined E_F is read from trace file).

Optional arguments:
    --condtens FILE   Path to *.condtens file  [default: PREFIX/interpolation.condtens]
    --trace    FILE   Path to *.trace file. If given, refined E_F (N=0 row) is used
                      instead of E_FERMI. Units in trace are Ry → converted automatically.
    --mu MIN MAX      Chemical-potential window in eV  [default: -1.5 1.5]
    --temp T [T …]    Temperatures to plot in K  [default: 300 500 700 900]
    --tau TAU         Relaxation time in seconds. If omitted, plots τ-normalised.
    --tau-fs TAU_FS   Same as --tau but in femtoseconds (e.g. --tau-fs 10 = 10 fs).
    --rescale LZ DEFF Apply 2D thickness rescaling: Lz and d_eff both in Angstrom.

Examples:
    python plot_btp2.py ZrClBr 5.2989
    python plot_btp2.py ZrClBr 5.2989 --tau-fs 10
    python plot_btp2.py ZrClBr 5.2989 --trace ./ZrClBr/interpolation.trace
    python plot_btp2.py ZrFBr  4.9365 --mu -2.0 2.0 --temp 300 700
    python plot_btp2.py ZrFCl  4.6109 --tau 1e-14 --rescale 20.0 6.15
    python plot_btp2.py ZrClBr 5.2989 --condtens ./ZrClBr/interpolation.condtens
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

RY_TO_EV = 13.605693122994


# ── helpers ────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot BoltzTraP2 thermoelectric properties.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("PREFIX",   help="Material prefix (e.g. ZrClBr)")
    p.add_argument("E_FERMI",  type=float, nargs="?", default=None,
                   help="Fermi/midgap energy in eV  [(VBM+CBM)/2]. If omitted, extracted from condtens.")

    p.add_argument("--condtens", metavar="FILE", default=None,
                   help="Path to *.condtens file  "
                        "[default: PREFIX/interpolation.condtens]")
    p.add_argument("--mu", nargs=2, type=float, metavar=("MIN", "MAX"),
                   default=[-1.5, 1.5],
                   help="Chemical-potential window in eV  [default: -1.5 1.5]")
    p.add_argument("--temp", nargs="+", type=float,
                   metavar="T", default=[300.0, 500.0, 700.0, 900.0],
                   help="Temperatures to plot in K  [default: 300 500 700 900]")
    p.add_argument("--tau", type=float, default=None, metavar="TAU",
                   help="Relaxation time in seconds  "
                        "[default: None → plot τ-normalised]")
    p.add_argument("--tau-fs", type=float, default=None, metavar="TAU_FS",
                   help="Relaxation time in femtoseconds  (e.g. --tau-fs 10 = 10 fs = 1e-14 s)")
    p.add_argument("--trace", metavar="FILE", default=None,
                   help="Path to *.trace file. Refined E_F (N≈0 row, in Ry) overrides E_FERMI.")
    p.add_argument("--rescale", nargs=2, type=float,
                   metavar=("LZ", "DEFF"), default=None,
                   help="2D thickness rescaling: Lz and d_eff in Angstrom")

    return p.parse_args()


def read_condtens(filename):
    """
    Read a BoltzTraP2 *.condtens file.

    Columns used here:
      0   Ef[Ry]
      1   T[K]
      2   N[e/uc]
      3..11   sigma/tau tensor   (xx yx zx  xy yy zy  xz yz zz)
      12..20  Seebeck tensor
      21..29  kappa_e/tau tensor
    """
    data = np.loadtxt(filename, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 30:
        raise RuntimeError(
            f"{filename} has {data.shape[1]} columns; at least 30 expected."
        )
    return data


def in_plane_average(data, col_xx, col_yy):
    return 0.5 * (data[:, col_xx] + data[:, col_yy])


def extract_efermi_from_data(data, source_name="trace"):
    if data.ndim == 1:
        data = data.reshape(1, -1)
    # Find row with smallest |N| at the first available temperature
    first_T   = data[0, 1]
    mask_T    = np.isclose(data[:, 1], first_T)
    sub       = data[mask_T]
    idx_min   = np.argmin(np.abs(sub[:, 2]))   # column 2 = N[e/uc]
    ef_ry     = sub[idx_min, 0]
    ef_ev     = ef_ry * RY_TO_EV
    print(f"[{source_name}] Refined E_F = {ef_ry:.8f} Ry = {ef_ev:.6f} eV  (N = {sub[idx_min,2]:.4e} e/uc)")
    return ef_ev


def read_refined_efermi(trace_file):
    """
    Read the refined Fermi energy from a BoltzTraP2 *.trace file.
    """
    data = np.loadtxt(trace_file, comments="#")
    return extract_efermi_from_data(data, source_name="trace")


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    PREFIX          = args.PREFIX
    FERMI_ENERGY_EV = args.E_FERMI
    MU_WINDOW_EV    = tuple(args.mu)
    TEMPERATURES_K  = args.temp

    # --tau-fs takes priority; convert fs → s
    if args.tau_fs is not None and args.tau is not None:
        sys.exit("ERROR: use either --tau or --tau-fs, not both.")
    if args.tau_fs is not None:
        TAU_SECONDS = args.tau_fs * 1.0e-15
        print(f"tau = {args.tau_fs} fs = {TAU_SECONDS:.3e} s")
    else:
        TAU_SECONDS = args.tau

    # Condtens path: explicit → per-material subdir → current dir
    if args.condtens:
        CONDTENS_FILE = Path(args.condtens)
    elif Path(f"{PREFIX}/interpolation.condtens").exists():
        CONDTENS_FILE = Path(f"{PREFIX}/interpolation.condtens")
    else:
        CONDTENS_FILE = Path("./interpolation.condtens")

    if not CONDTENS_FILE.exists():
        sys.exit(
            f"ERROR: {CONDTENS_FILE} not found.\n"
            "Run 'btp2 integrate' first, or pass --condtens <path>."
        )

    # ── load data ──────────────────────────────────────────────────────────────
    data        = read_condtens(CONDTENS_FILE)

    # --trace: override E_FERMI with refined value from trace file (units: Ry → eV)
    if args.trace:
        trace_path = Path(args.trace)
        if not trace_path.exists():
            sys.exit(f"ERROR: trace file {trace_path} not found.")
        FERMI_ENERGY_EV = read_refined_efermi(trace_path)
    elif FERMI_ENERGY_EV is not None:
        print(f"E_F = {FERMI_ENERGY_EV:.6f} eV  (from command-line argument)")
    else:
        FERMI_ENERGY_EV = extract_efermi_from_data(data, source_name="condtens")

    # 2D rescaling
    thickness_factor = 1.0
    rescale_label    = "raw supercell normalization"
    if args.rescale:
        lz, deff         = args.rescale
        thickness_factor = lz / deff
        rescale_label    = f"rescaled Lz/d_eff = {lz:g}/{deff:g}"

    OUTPUT_PNG = Path(f"{PREFIX}_thermoelectric_vs_mu.png")

    mu_rel_ev   = data[:, 0] * RY_TO_EV - FERMI_ENERGY_EV
    temperature = data[:, 1]

    sigma_over_tau = in_plane_average(data, 3, 7)  * thickness_factor
    seebeck_v_k    = in_plane_average(data, 12, 16)
    kappa_over_tau = in_plane_average(data, 21, 25) * thickness_factor

    seebeck_uv_k   = seebeck_v_k * 1.0e6

    # ── apply tau if given ─────────────────────────────────────────────────────
    if TAU_SECONDS is None:
        sigma_plot = sigma_over_tau
        kappa_plot = kappa_over_tau
        sigma_ylabel = r"$\sigma_{\parallel}/\tau$ [$(\Omega\,m\,s)^{-1}$]"
        kappa_ylabel = r"$\kappa_{e,\parallel}/\tau$ [W m$^{-1}$ K$^{-1}$ s$^{-1}$]"
        tau_label    = "per relaxation time"
    else:
        sigma_plot = sigma_over_tau * TAU_SECONDS
        kappa_plot = kappa_over_tau * TAU_SECONDS
        sigma_ylabel = r"$\sigma_{\parallel}$ [S m$^{-1}$]"
        kappa_ylabel = r"$\kappa_{e,\parallel}$ [W m$^{-1}$ K$^{-1}$]"
        tau_label    = f"tau = {TAU_SECONDS:.3e} s"

    # ── plot ───────────────────────────────────────────────────────────────────
    print(f"Temperatures available in file:")
    print(np.unique(temperature))

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0), sharex=True)
    ax_s, ax_sigma, ax_kappa = axes.ravel()

    plotted_any = False
    for target_T in TEMPERATURES_K:
        mask = np.isclose(temperature, target_T, rtol=0.0, atol=1.0e-8)
        if not np.any(mask):
            print(f"WARNING: T = {target_T:g} K not found in file. Skipping.")
            continue

        x   = mu_rel_ev[mask]
        s   = seebeck_uv_k[mask]
        sig = sigma_plot[mask]
        kap = kappa_plot[mask]

        order  = np.argsort(x)
        x, s, sig, kap = x[order], s[order], sig[order], kap[order]

        window = (x >= MU_WINDOW_EV[0]) & (x <= MU_WINDOW_EV[1])
        if not np.any(window):
            continue
        x, s, sig, kap = x[window], s[window], sig[window], kap[window]

        label = f"{target_T:g} K"
        ax_s.plot(x, s, label=label)
        ax_sigma.plot(x, sig, label=label)
        ax_kappa.plot(x, kap, label=label)
        plotted_any = True

    if not plotted_any:
        sys.exit("ERROR: No curves plotted. Check --temperatures and --mu arguments.")

    for ax in axes.ravel():
        ax.axvline(0.0, linewidth=0.9, linestyle="--", color="k")
        ax.set_xlim(*MU_WINDOW_EV)
        ax.legend(frameon=False)

    ax_s.axhline(0.0, linewidth=0.8, linestyle=":", color="k")
    ax_s.set_ylabel(r"$S_{\parallel}$ [$\mu$V/K]")
    ax_s.set_title("Seebeck coefficient")
    ax_s.set_xlabel(r"$\mu - E_F$ (eV)")

    ax_sigma.set_ylabel(sigma_ylabel)
    ax_sigma.set_title("Electrical conductivity")
    ax_sigma.set_xlabel(r"$\mu - E_F$ (eV)")

    ax_kappa.set_ylabel(kappa_ylabel)
    ax_kappa.set_title("Electronic thermal conductivity")
    ax_kappa.set_xlabel(r"$\mu - E_F$ (eV)")

    fig.suptitle(
        rf"{PREFIX} — Thermoelectric Properties  [{tau_label}]"
    )
    fig.tight_layout()

    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")

    print()
    print(f"Prefix          : {PREFIX}")
    print(f"Reference E_F   : {FERMI_ENERGY_EV:.6f} eV")
    print(f"condtens file   : {CONDTENS_FILE}")
    print(f"Normalization   : {rescale_label}")
    print(f"PNG written to  : {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
