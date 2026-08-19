#!/usr/bin/env python3
"""Plot Quantum ESPRESSO DOS and orbital-resolved PDOS with SciencePlots.

The script reads ``<prefix>.dos`` and the files produced by ``projwfc.x``:
``<prefix>.pdos.pdos_atm#*(*)_wfc#*(*)``.  PDOS contributions are summed by
element and angular momentum, so equivalent atoms are shown as a single curve.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scienceplots  # noqa: F401  (registers the SciencePlots styles)


PDOS_PATTERN = re.compile(
    r"\.pdos_atm#(?P<atom>\d+)\((?P<element>[^)]+)\)_wfc#\d+\((?P<orbital>[spdf])\)$"
)
ORBITAL_ORDER = {"s": 0, "p": 1, "d": 2, "f": 3}
COLORS = {
    # Paul Tol-inspired, high-contrast and colour-vision-friendly palette.
    ("Zr", "s"): "#4477AA",  # blue
    ("Zr", "p"): "#44AA99",  # teal
    ("Zr", "d"): "#EE6677",  # coral red
    ("Cl", "s"): "#228833",  # green
    ("Cl", "p"): "#AA3377",  # purple
    ("Br", "s"): "#EE7733",  # orange
    ("Br", "p"): "#8C564B",  # brown
}
LINESTYLES = {
    ("Zr", "s"): "-",
    ("Zr", "p"): (0, (6, 2)),
    ("Zr", "d"): "-",
    ("Cl", "s"): (0, (4, 1.5, 1, 1.5)),
    ("Cl", "p"): "-",
    ("Br", "s"): (0, (2, 1.4)),
    ("Br", "p"): (0, (7, 2, 1.5, 2)),
}


def str2bool(value: str | bool) -> bool:
    """Match the boolean CLI convention used by the original plotdos.py."""
    if isinstance(value, bool):
        return value
    if value.lower() in {"yes", "true", "t", "y", "1"}:
        return True
    if value.lower() in {"no", "false", "f", "n", "0"}:
        return False
    raise argparse.ArgumentTypeError("Boolean value expected.")


def read_table(path: Path) -> np.ndarray:
    """Read a whitespace-separated QE table while ignoring its comment header."""
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Format data tidak dikenali: {path}")
    return data


def fermi_from_dos_header(path: Path) -> float:
    """Read the Fermi/reference energy written by dos.x in the first line."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        header = handle.readline()
    match = re.search(r"EFermi\s*=\s*([-+0-9.eEdD]+)", header, re.IGNORECASE)
    if not match:
        raise ValueError(f"EFermi tidak ditemukan pada header {path}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))


def collect_pdos(data_dir: Path, prefix: str) -> tuple[np.ndarray, dict[tuple[str, str], np.ndarray]]:
    """Sum each file's l-resolved DOS (column 2) by element and orbital."""
    grouped: dict[tuple[str, str], np.ndarray] = defaultdict(lambda: None)
    energy: np.ndarray | None = None
    paths = sorted(data_dir.glob(f"{prefix}.pdos.pdos_atm#*_wfc#*"))
    if not paths:
        raise FileNotFoundError(f"File PDOS tidak ditemukan di {data_dir}")

    for path in paths:
        match = PDOS_PATTERN.search(path.name)
        if not match:
            continue
        table = read_table(path)
        if energy is None:
            energy = table[:, 0]
        elif table.shape[0] != energy.size or not np.allclose(table[:, 0], energy):
            raise ValueError(f"Grid energi PDOS tidak konsisten: {path}")

        key = (match.group("element"), match.group("orbital"))
        if grouped[key] is None:
            grouped[key] = np.zeros_like(table[:, 1])
        grouped[key] += table[:, 1]

    if energy is None or not grouped:
        raise ValueError("Tidak ada file PDOS atom/orbital yang dapat dibaca")
    return energy, dict(grouped)


def plot_dos_pdos(
    prefix: str,
    output: str | Path | None = None,
    fermi_in: float | None = None,
    data_file: str | Path | None = None,
    pdos_dir: str | Path | None = None,
    Emin: float = -8.0,
    Estep: float = 2.0,
    Emax: float = 8.0,
    Ymax: float = 20.0,
    dpi: int = 300,
    dark: bool = False,
    no_fill: bool = False,
    shift: bool = True,
) -> None:
    # Variable names and defaults intentionally follow the original plotdos.py.
    dos_path = Path(data_file) if data_file else Path(f"{prefix}.dos")
    if not dos_path.exists():
        raise FileNotFoundError(f"File DOS tidak ditemukan: {dos_path}")

    pdos_path = Path(pdos_dir) if pdos_dir else dos_path.resolve().parent
    out_file = Path(output) if output else Path(f"{prefix}_dos_pdos.png")
    separate_file = out_file.with_name(f"{out_file.stem}_separate{out_file.suffix}")
    stacked_file = out_file.with_name(f"{out_file.stem}_stacked{out_file.suffix}")

    dos = read_table(dos_path)
    pdos_energy, grouped = collect_pdos(pdos_path, prefix)
    reference = fermi_from_dos_header(dos_path) if fermi_in is None else fermi_in
    shift_value = reference if shift else 0.0
    dos_energy = dos[:, 0] - shift_value
    pdos_energy = pdos_energy - shift_value

    style = ["science", "no-latex"]
    if dark:
        style.append("dark_background")
    total_color = "white" if dark else "black"
    reference_color = "0.8" if dark else "0.25"
    legend_face = "0.08" if dark else "white"

    with plt.style.context(style):
        fig, ax = plt.subplots(figsize=(7.0, 5.2))

        ax.plot(
            dos_energy,
            dos[:, 1],
            color=total_color,
            lw=2.0,
            alpha=0.68,
            label="Total DOS",
            zorder=2,
        )
        if not no_fill:
            ax.fill_between(dos_energy, dos[:, 1], color="0.70", alpha=0.12, zorder=0)

        keys = sorted(
            grouped,
            key=lambda key: (key[0] != "Zr", key[0], ORBITAL_ORDER.get(key[1], 99)),
        )
        fallback_colors = plt.colormaps["tab10"]
        for index, key in enumerate(keys):
            element, orbital = key
            color = COLORS.get(key, fallback_colors(index % 10))
            ax.plot(
                pdos_energy,
                grouped[key],
                lw=2.0,
                color=color,
                linestyle=LINESTYLES.get(key, "-"),
                alpha=0.98,
                label=rf"{element} ${orbital}$",
                zorder=3,
            )

        # Confirm that the requested energy window contains both datasets.
        dos_mask = (dos_energy >= Emin) & (dos_energy <= Emax)
        pdos_mask = (pdos_energy >= Emin) & (pdos_energy <= Emax)
        if not dos_mask.any() or not pdos_mask.any():
            raise ValueError("Tidak ada data di dalam rentang energi yang dipilih")

        ax.axvline(
            reference - shift_value,
            color=reference_color,
            lw=1.4,
            ls="--",
            zorder=1,
        )
        ax.set_xlim(Emin, Emax)
        ax.set_ylim(0, Ymax)
        ax.set_ylabel(r"DOS (states eV$^{-1}$)", fontsize=17)
        ax.set_xticks(
            np.arange(
                np.ceil(Emin / Estep) * Estep,
                np.floor(Emax / Estep) * Estep + Estep * 0.01,
                Estep,
            )
        )
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.tick_params(axis="both", which="major", labelsize=14, width=1.2, length=6)
        ax.tick_params(axis="both", which="minor", width=1.0, length=3.5)

        ax.legend(
            frameon=True,
            facecolor=legend_face,
            edgecolor="none",
            framealpha=0.96,
            ncol=4 if len(keys) + 1 > 6 else 3,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.01),
            fontsize=12.5,
            columnspacing=1.2,
            handlelength=2.1,
        )

        title_prefix = prefix.replace("Cl2", r"Cl$_2$")
        ax.set_title(rf"DOS and PDOS of {title_prefix}", fontsize=21, pad=55)
        energy_label = r"$E-E_{\mathrm{VBM}}$ (eV)" if shift else r"Energy (eV)"
        ax.set_xlabel(energy_label, fontsize=17)

        fig.tight_layout()
        out_file.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_file, dpi=dpi, bbox_inches="tight")
        if out_file.suffix.lower() != ".pdf":
            fig.savefig(out_file.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)

        # A second output with separate DOS and PDOS panels.
        fig_sep, (ax_dos, ax_pdos) = plt.subplots(2, 1, figsize=(7.0, 6.8), sharex=True)
        ax_dos.plot(
            dos_energy,
            dos[:, 1],
            color=total_color,
            lw=2.0,
            alpha=0.78,
            label="Total DOS",
            zorder=2,
        )
        if not no_fill:
            ax_dos.fill_between(dos_energy, dos[:, 1], color="0.70", alpha=0.12, zorder=0)

        for index, key in enumerate(keys):
            element, orbital = key
            color = COLORS.get(key, fallback_colors(index % 10))
            ax_pdos.plot(
                pdos_energy,
                grouped[key],
                lw=2.0,
                color=color,
                linestyle=LINESTYLES.get(key, "-"),
                alpha=0.98,
                label=rf"{element} ${orbital}$",
                zorder=3,
            )

        for panel in (ax_dos, ax_pdos):
            panel.axvline(
                reference - shift_value,
                color=reference_color,
                lw=1.4,
                ls="--",
                zorder=1,
            )
            panel.set_xlim(Emin, Emax)
            panel.set_ylim(0, Ymax)
            panel.set_ylabel(r"DOS (states eV$^{-1}$)", fontsize=17)
            panel.set_xticks(
                np.arange(
                    np.ceil(Emin / Estep) * Estep,
                    np.floor(Emax / Estep) * Estep + Estep * 0.01,
                    Estep,
                )
            )
            panel.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
            panel.tick_params(axis="both", which="major", labelsize=14, width=1.2, length=6)
            panel.tick_params(axis="both", which="minor", width=1.0, length=3.5)

        ax_dos.legend(frameon=False, loc="upper right", fontsize=14)
        ax_pdos.legend(
            frameon=True,
            facecolor=legend_face,
            edgecolor="none",
            framealpha=0.96,
            ncol=3 if len(keys) > 5 else 2,
            loc="upper right",
            fontsize=13,
            columnspacing=1.4,
            handlelength=2.4,
        )
        ax_dos.set_title(rf"DOS and PDOS of {title_prefix}", fontsize=21, pad=12)
        ax_pdos.set_title("Orbital-resolved PDOS", loc="left", fontsize=16, pad=7)
        ax_pdos.set_xlabel(energy_label, fontsize=17)
        fig_sep.align_ylabels((ax_dos, ax_pdos))
        fig_sep.tight_layout()
        fig_sep.savefig(separate_file, dpi=dpi, bbox_inches="tight")
        if separate_file.suffix.lower() != ".pdf":
            fig_sep.savefig(separate_file.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig_sep)

        # A third, qualitative mirrored-panel view inspired by stacked DOS
        # figures. Each strip is normalized independently so weak orbital
        # channels remain visible; the original absolute-scale plots above are
        # retained for quantitative comparison.
        stacked_series = [("Total", dos_energy, dos[:, 1], "0.45", "-")]
        for index, key in enumerate(keys):
            element, orbital = key
            stacked_series.append(
                (
                    rf"{element} ${orbital}$",
                    pdos_energy,
                    grouped[key],
                    COLORS.get(key, fallback_colors(index % 10)),
                    LINESTYLES.get(key, "-"),
                )
            )

        n_panels = len(stacked_series)
        fig_stack, stack_axes = plt.subplots(
            n_panels,
            1,
            figsize=(7.0, 1.05 * n_panels + 1.25),
            sharex=True,
            gridspec_kw={"hspace": 0.0},
        )
        stack_axes = np.atleast_1d(stack_axes)
        for panel, (label, energy, values, color, linestyle) in zip(stack_axes, stacked_series):
            mask = (energy >= Emin) & (energy <= Emax)
            peak = float(np.max(np.abs(values[mask])))
            normalized = values / peak if peak > 0 else np.zeros_like(values)

            panel.fill_between(energy, normalized, -normalized, color=color, alpha=0.92)
            panel.plot(energy, normalized, color=color, lw=1.0, linestyle=linestyle)
            panel.plot(energy, -normalized, color=color, lw=1.0, linestyle=linestyle)
            panel.axhline(0.0, color="0.55", lw=0.7, zorder=0)
            panel.axvline(
                reference - shift_value,
                color="0.55",
                lw=0.9,
                ls=":",
                zorder=4,
            )
            panel.set_xlim(Emin, Emax)
            panel.set_ylim(-1.08, 1.08)
            panel.set_yticks([])
            panel.text(
                0.985,
                0.77,
                label,
                transform=panel.transAxes,
                ha="right",
                va="center",
                fontsize=14,
                bbox={"facecolor": legend_face, "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
            )
            panel.tick_params(
                axis="x",
                which="major",
                direction="in",
                top=True,
                bottom=True,
                length=5,
                width=1.0,
                labelsize=14,
            )
            panel.tick_params(
                axis="x",
                which="minor",
                direction="in",
                top=True,
                bottom=True,
                length=2.5,
                width=0.8,
            )
            panel.xaxis.set_minor_locator(ticker.AutoMinorLocator(4))

        for panel in stack_axes[:-1]:
            panel.tick_params(labelbottom=False)
        stack_axes[-1].set_xticks(
            np.arange(
                np.ceil(Emin / Estep) * Estep,
                np.floor(Emax / Estep) * Estep + Estep * 0.01,
                Estep,
            )
        )
        stack_axes[-1].set_xlabel(energy_label, fontsize=17)
        fig_stack.supylabel("Normalized DOS", fontsize=17, x=0.015)
        fig_stack.suptitle(title_prefix, fontsize=21, y=0.995)
        fig_stack.subplots_adjust(left=0.14, right=0.98, bottom=0.09, top=0.94)
        fig_stack.savefig(stacked_file, dpi=dpi, bbox_inches="tight")
        if stacked_file.suffix.lower() != ".pdf":
            fig_stack.savefig(stacked_file.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig_stack)

    print(f"Reference energy : {reference:.3f} eV")
    print(f"Applied shift    : {-shift_value:+.3f} eV")
    print(f"Energy window    : {Emin:.2f} to {Emax:.2f} eV (step {Estep:g} eV)")
    print(f"Y-axis range     : 0 to {Ymax:g} states/eV")
    print(f"PDOS channels    : {', '.join(f'{e}-{o}' for e, o in keys)}")
    print(f"Saved            : {out_file.resolve()}")
    if out_file.suffix.lower() != ".pdf":
        print(f"Saved            : {out_file.with_suffix('.pdf').resolve()}")
    print(f"Saved            : {separate_file.resolve()}")
    if separate_file.suffix.lower() != ".pdf":
        print(f"Saved            : {separate_file.with_suffix('.pdf').resolve()}")
    print(f"Saved            : {stacked_file.resolve()}")
    if stacked_file.suffix.lower() != ".pdf":
        print(f"Saved            : {stacked_file.with_suffix('.pdf').resolve()}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plot_dos_pdos.py",
        description="Plot Quantum ESPRESSO DOS and orbital-resolved PDOS with SciencePlots",
    )
    parser.add_argument("prefix", help="QE prefix")
    parser.add_argument("-o", "--output", metavar="FILE", help="Default: <prefix>_dos_pdos.png")
    parser.add_argument("-f", "--fermi", metavar="EV", type=float, help="Override Fermi/VBM energy")
    parser.add_argument("--data", metavar="PATH", help="DOS data file (default: <prefix>.dos)")
    parser.add_argument(
        "--pdos-dir",
        "--data-dir",
        dest="pdos_dir",
        metavar="DIR",
        help="PDOS folder (default: folder containing the DOS file)",
    )
    parser.add_argument("--dpi", metavar="INT", type=int, default=300, help="PNG DPI (default: 300)")
    parser.add_argument("--ymax", metavar="VALUE", type=float, default=20.0, help="Y-axis maximum (default: 20)")
    parser.add_argument("--dark", action="store_true", help="Use a dark background")
    parser.add_argument("--no-fill", action="store_true", help="Do not shade the total DOS")
    parser.add_argument(
        "--shift",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Shift the reference energy to 0 eV (default: True)",
    )
    parser.add_argument("Emin", nargs="?", type=float, default=-8.0)
    parser.add_argument("Estep", nargs="?", type=float, default=2.0)
    parser.add_argument("Emax", nargs="?", type=float, default=8.0)
    return parser


def main() -> None:
    # Allow options and the optional Emin/Estep/Emax values to be interleaved.
    args = build_parser().parse_intermixed_args()
    if args.Emin >= args.Emax:
        raise SystemExit("Emin harus lebih kecil daripada Emax")
    if args.Estep <= 0:
        raise SystemExit("Estep harus lebih besar daripada 0")
    if args.ymax <= 0:
        raise SystemExit("--ymax harus lebih besar daripada 0")
    plot_dos_pdos(
        prefix=args.prefix,
        output=args.output,
        fermi_in=args.fermi,
        data_file=args.data,
        pdos_dir=args.pdos_dir,
        Emin=args.Emin,
        Estep=args.Estep,
        Emax=args.Emax,
        Ymax=args.ymax,
        dpi=args.dpi,
        dark=args.dark,
        no_fill=args.no_fill,
        shift=args.shift,
    )


if __name__ == "__main__":
    main()
