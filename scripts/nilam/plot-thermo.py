#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
CONDTENS_FILE = SCRIPT_DIR / "interpolation.condtens"

# Intrinsic reference energy from the SCF band edges, in eV:
#
# E_F = (E_VBM + E_CBM) / 2

# SCF
VBM_ENERGY_EV = 3.7915
CBM_ENERGY_EV = 4.8247
FERMI_ENERGY_EV = 0.5 * (VBM_ENERGY_EV + CBM_ENERGY_EV)

# NSCF
# VBM_ENERGY_EV = 3.866500
# CBM_ENERGY_EV = 4.824800
# FERMI_ENERGY_EV = 4.824800

TEMPERATURES_K = [300.0, 600.0, 900.0]
MU_WINDOW_EV = (-1.5, 1.5)
PLOT_POWER_FACTOR = False

# Leave this as None to plot sigma/tau, kappa_e/tau, and PF/tau.
# Enter a value only if the relaxation time has a defensible basis.
TAU_SECONDS = 10e-15


# Optional 2D thickness normalization.
APPLY_2D_RESCALE = False
LZ_ANGSTROM = 20.0
D_EFFECTIVE_ANGSTROM = 6.15

OUTPUT_PNG = SCRIPT_DIR / "figures/ZrCl2_thermoelectric_vs_mu.png"
OUTPUT_PDF = SCRIPT_DIR / "figures/ZrCl2_thermoelectric_vs_mu.pdf"

RY_TO_EV = 13.605693122994


def read_condtens(filename):
    """
    Read a BoltzTraP2 *.condtens file.

    Columns used here:
      0   Ef[Ry]
      1   T[K]
      2   N[e/uc]

      3..11   sigma/tau tensor
      12..20  Seebeck tensor
      21..29  electronic thermal conductivity/tau tensor

    Tensor order is:
      xx, yx, zx, xy, yy, zy, xz, yz, zz
    """
    data = np.loadtxt(filename, comments="#")

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 30:
        raise RuntimeError(
            f"{filename} has {data.shape[1]} columns. "
            "At least 30 columns are expected."
        )

    return data


def in_plane_average(data, col_xx, col_yy):
    return 0.5 * (data[:, col_xx] + data[:, col_yy])


def main():
    if not CONDTENS_FILE.exists():
        raise FileNotFoundError(
            f"{CONDTENS_FILE} was not found. Run 'btp2 integrate' first."
        )

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    data = read_condtens(CONDTENS_FILE)

    mu_abs_ev = data[:, 0] * RY_TO_EV
    mu_rel_ev = mu_abs_ev - FERMI_ENERGY_EV
    temperature = data[:, 1]

    sigma_xx_over_tau = data[:, 3]
    sigma_yy_over_tau = data[:, 7]
    seebeck_xx_v_per_k = data[:, 12]
    seebeck_yy_v_per_k = data[:, 16]

    sigma_over_tau = 0.5 * (sigma_xx_over_tau + sigma_yy_over_tau)
    seebeck_v_per_k = 0.5 * (
        seebeck_xx_v_per_k + seebeck_yy_v_per_k
    )
    kappa_over_tau = in_plane_average(data, 21, 25)

    thickness_factor = 1.0

    if APPLY_2D_RESCALE:
        thickness_factor = LZ_ANGSTROM / D_EFFECTIVE_ANGSTROM
        sigma_xx_over_tau = sigma_xx_over_tau * thickness_factor
        sigma_yy_over_tau = sigma_yy_over_tau * thickness_factor
        sigma_over_tau = sigma_over_tau * thickness_factor
        kappa_over_tau = kappa_over_tau * thickness_factor

    # Power factor untuk setiap arah dihitung dengan:
    #
    #     PF_ii / tau = S_ii**2 * (sigma_ii / tau)
    #
    # S harus tetap dalam V/K pada perhitungan ini, bukan microvolt/K.
    # Dengan satuan dari file condtens, hasilnya mempunyai satuan
    # W m^-1 K^-2 s^-1. Nilai PF sebenarnya baru diperoleh setelah
    # dikalikan dengan relaxation time, tau.
    #
    # Hitung PF_xx dan PF_yy terlebih dahulu, kemudian ambil rata-rata
    # bidang. Secara umum mean(S**2 * sigma) tidak sama dengan
    # mean(S)**2 * mean(sigma).
    pf_over_tau = 0.5 * (
        seebeck_xx_v_per_k**2 * sigma_xx_over_tau
        + seebeck_yy_v_per_k**2 * sigma_yy_over_tau
    )
    seebeck_uv_per_k = seebeck_v_per_k * 1.0e6

    if TAU_SECONDS is None:
        sigma_plot = sigma_over_tau
        kappa_plot = kappa_over_tau
        pf_plot = pf_over_tau

        sigma_ylabel = r"$\sigma_{\parallel}/\tau$ [$(\Omega\,m\,s)^{-1}$]"
        kappa_ylabel = (
            r"$\kappa_{e,\parallel}/\tau$ "
            r"[W m$^{-1}$ K$^{-1}$ s$^{-1}$]"
        )
        pf_ylabel = (
            r"$PF_{\parallel}/\tau$ "
            r"[W m$^{-1}$ K$^{-2}$ s$^{-1}$]"
        )
    else:
        sigma_plot = sigma_over_tau * TAU_SECONDS
        kappa_plot = kappa_over_tau * TAU_SECONDS
        pf_plot = pf_over_tau * TAU_SECONDS

        sigma_ylabel = r"$\sigma_{\parallel}$ [S m$^{-1}$]"
        kappa_ylabel = r"$\kappa_{e,\parallel}$ [W m$^{-1}$ K$^{-1}$]"
        pf_ylabel = r"$PF_{\parallel}$ [W m$^{-1}$ K$^{-2}$]"

    if PLOT_POWER_FACTOR:
        fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), sharex=True)
        ax_s, ax_sigma, ax_kappa, ax_pf = axes.ravel()
    else:
        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), sharex=True)
        ax_s, ax_sigma, ax_kappa = axes.ravel()
        ax_pf = None

    plotted_any = False

    print("Temperatures available in the file:")
    print(np.unique(temperature))

    for target_T in TEMPERATURES_K:
        mask = np.isclose(
            temperature,
            target_T,
            rtol=0.0,
            atol=1.0e-8,
        )

        if not np.any(mask):
            print(f"WARNING: T = {target_T:g} K was not found. Skipping.")
            continue

        x = mu_rel_ev[mask]
        s = seebeck_uv_per_k[mask]
        sig = sigma_plot[mask]
        kap = kappa_plot[mask]
        pf = pf_plot[mask]

        order = np.argsort(x)

        x = x[order]
        s = s[order]
        sig = sig[order]
        kap = kap[order]
        pf = pf[order]

        window = (
            (x >= MU_WINDOW_EV[0])
            &
            (x <= MU_WINDOW_EV[1])
        )

        if not np.any(window):
            continue

        x = x[window]
        s = s[window]
        sig = sig[window]
        kap = kap[window]
        pf = pf[window]

        label = f"{target_T:g} K"

        ax_s.plot(x, s, label=label)
        ax_sigma.plot(x, sig, label=label)
        ax_kappa.plot(x, kap, label=label)
        # Panel power factor dinonaktifkan melalui PLOT_POWER_FACTOR.
        if PLOT_POWER_FACTOR:
            ax_pf.plot(x, pf, label=label)

        plotted_any = True

    if not plotted_any:
        raise RuntimeError(
            "No curves were plotted. Check TEMPERATURES_K "
            "and MU_WINDOW_EV."
        )

    for ax in axes.ravel():
        ax.axvline(0.0, linewidth=0.9, linestyle="--")
        ax.set_xlim(*MU_WINDOW_EV)
        ax.grid(False)

    ax_s.axhline(0.0, linewidth=0.8, linestyle=":")
    ax_s.set_ylabel(r"$S_{\parallel}$ [$\mu$V/K]")
    ax_s.set_title("Seebeck coefficient")

    ax_sigma.set_ylabel(sigma_ylabel)
    ax_sigma.set_title("Electrical conductivity")

    ax_kappa.set_ylabel(kappa_ylabel)
    ax_kappa.set_title("Electronic thermal conductivity")

    if PLOT_POWER_FACTOR:
        ax_kappa.set_xlabel(r"$\mu-E_F$ (eV)")

        # Plot power factor bidang terhadap chemical potential relatif
        # terhadap energi Fermi intrinsik. Karena TAU_SECONDS masih None,
        # kurva yang ditampilkan adalah PF/tau, bukan PF absolut.
        ax_pf.set_ylabel(pf_ylabel)
        ax_pf.set_title("Power factor")
        ax_pf.set_xlabel(r"$\mu-E_F$ (eV)")
    else:
        for ax in axes.ravel():
            ax.set_xlabel(r"$\mu-E_F$ (eV)")

    for ax in axes.ravel():
        ax.legend(frameon=False)

    normalization_text = "raw supercell normalization"

    if APPLY_2D_RESCALE:
        normalization_text = (
            f"rescaled with Lz/d_eff = "
            f"{LZ_ANGSTROM:g}/{D_EFFECTIVE_ANGSTROM:g}"
        )

    tau_text = (
        "per relaxation time"
        if TAU_SECONDS is None
        else f"tau = {TAU_SECONDS:.3e} s"
    )

    fig.suptitle(
        f"ZrCl$_2$ — {tau_text}, {normalization_text}"
    )
    fig.tight_layout()

    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")

    print()
    print(f"Reference E_F = {FERMI_ENERGY_EV:.8f} eV")
    print(f"Thickness normalization factor = {thickness_factor:.8f}")
    print(f"PNG figure written to: {OUTPUT_PNG}")
    print(f"PDF figure written to: {OUTPUT_PDF}")

    plt.show()


if __name__ == "__main__":
    main()
