# Electronic, Elastic, and Thermoelectric Transport Properties of 2D $\mathrm{Zr}X_2$ and Janus $\mathrm{Zr}XY$ Monolayers

This repository hosts collaborative first-principles calculations (Quantum ESPRESSO, BoltzTraP2, and Phonopy/Phono3py) investigating the electronic structure, acoustic phonon scattering (Deformation Potential Theory), and thermoelectric transport properties of 2D transition metal dihalides ($\mathrm{Zr}X_2: \mathrm{ZrF}_2, \mathrm{ZrCl}_2, \mathrm{ZrBr}_2$) and Janus monolayers ($\mathrm{Zr}XY: \mathrm{ZrFCl}, \mathrm{ZrFBr}, \mathrm{ZrClBr}$).

---

## 📁 Repository Structure

```
ZrXY-draft/
├── README.md
├── .gitignore
│
├── ZrF2/                     <- Parent Difluoride Monolayer
│   ├── agna/
│   └── rizka/
│
├── ZrCl2/                    <- Parent Dichloride Monolayer
│   ├── nilam/
│   └── rizka/
│
├── ZrBr2/                    <- Parent Dibromide Monolayer
│   └── jauhar/
│
├── ZrClBr/                   <- Janus Monolayer (ZrClBr)
│   └── rizka/
│       ├── band/             <- Standard band structure (nscf, bands.x)
│       ├── bands2/           <- Super-dense k-point mesh & effective mass fitting
│       ├── btp2/             <- BoltzTraP2 transport interpolation & raw profiles
│       ├── dos/              <- Total and projected density of states (dos.x, projwfc.x)
│       ├── strain/           <- Uniaxial strain (dir_x, dir_y), C2D, & vacuum-aligned E1
│       └── phonon/           <- Phonon dispersion & harmonic stability (phonopy)
│
├── ZrFBr/                    <- Janus Monolayer (ZrFBr)
│   ├── jauhar/               <- Reference calculations (DPT benchmark, bands2, strain)
│   └── rizka/                <- Full production calculations (band, bands2, btp2, dos, pdos, strain, phonon)
│
├── ZrFCl/                    <- Janus Monolayer (ZrFCl)
│   └── rizka/                <- Full production calculations (band, bands2, btp2, dos, pdos, strain, phonon)
│
└── scripts/
    ├── nilam/                <- Helper conversion and plotting scripts
    └── rizka/                <- Master Slurm job runners, DPT solvers, and phonon setup
        ├── job.all_bands.sh            <- Batch runner for super-dense bands (all materials)
        ├── job.all_strains.sh          <- Batch runner for DPT strain + vacuum alignment
        ├── job.all_phonons_harmonic.sh <- Batch runner for harmonic phonon supercells
        ├── setup_harmonic_phonons.py   <- Automated phonopy-init displacement generator
        └── plot_pf_physical_dpt.py     <- Solver for physical Power Factor with DPT relaxation time
```

---

## Fundamental Electronic Properties & Fermi Energies

Fermi levels ($E_F$) defined at mid-gap $\frac{E_{\mathrm{VBM}} + E_{\mathrm{CBM}}}{2}$ from high-precision SCF/NSCF:

| Material | Type | Lattice $a\ (\text{\AA})$ | VBM $(\text{eV})$ | CBM $(\text{eV})$ | Band Gap $E_g\ (\text{eV})$ | Mid-Gap $E_F\ (\text{eV})$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **$\mathrm{ZrFCl}$** | Janus | $3.2589$ | $4.0173$ | $5.2918$ | **$1.2745$** (Indirect) | **$4.6546$** |
| **$\mathrm{ZrFBr}$** | Janus | $3.3444$ | $4.3671$ | $5.5612$ | **$1.1941$** (Indirect) | **$4.9642$** |
| **$\mathrm{ZrClBr}$**| Janus | $3.4891$ | $4.8960$ | $5.7991$ | **$0.9031$** (Indirect) | **$5.3476$** |
| **$\mathrm{ZrBr}_2$**| Parent| $3.3439$ | — | — | **$0.8120$** (Indirect) | — |

---

## Thermoelectric & Deformation Potential Theory (DPT) Summary

To overcome the Constant Relaxation Time Approximation (CRTA) in standard BoltzTraP2, intrinsic carrier relaxation times $\tau_0 \equiv \tau_{\mathrm{DPT}}$ are calculated using the 2D Bardeen–Shockley acoustic phonon scattering relation:

$$\tau(T) = \frac{2 \hbar^3 C_{2D}}{3 k_B T \, m^* \, (E_1)^2}$$

$$\mu(T) = \frac{e \, \tau(T)}{m^*} = \frac{2 e \hbar^3 C_{2D}}{3 k_B T \, (m^*)^2 \, (E_1)^2}$$

### Physical Parameters ($T = 300\ \text{K}$):

| Material | Carrier | $m^*/m_0$ | $C_{2D}\ (\text{N/m})$ | $\|E_1\|\ (\text{eV})$ | $\tau(300\text{ K})\ (\text{fs})$ | $\mu(300\text{ K})\ (\frac{\text{cm}^2}{\text{Vs}})$ | $PF_{\max}\ (\frac{\text{mW}}{\text{mK}^2})$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **$\mathrm{ZrFCl}$** | Electron ($n$) | $1.590$ | $100.03$ | $4.38$ | **$32.9$** | $36.4$ | **$34.84$** |
| **$\mathrm{ZrFCl}$** | Hole ($p$)     | $0.629$ | $100.03$ | $7.67$ | **$26.7$** | $74.7$ | **$13.46$** |
| **$\mathrm{ZrFBr}$** | Electron ($n$) | $1.245$ | $94.85$  | $4.87$ | **$30.6$** | $43.2$ | **$35.57$** |
| **$\mathrm{ZrFBr}$** | Hole ($p$)     | $0.707$ | $94.85$  | $7.56$ | **$19.4$** | $48.3$ | **$23.21$** |
| **$\mathrm{ZrClBr}$**| Electron ($n$) | $1.096$ | $78.16$  | $5.99$ | **$19.3$** | $31.0$ | **$23.04$** |
| **$\mathrm{ZrClBr}$**| Hole ($p$)     | $0.459$ | $78.16$  | $9.82$ | **$17.2$** | $65.8$ | **$20.78$** |

---

## Key Methodological Features

1. **Super-Dense Bands (`bands2/`)**:
   - Ultra-dense $k$-point sampling around the Valence Band Maximum (VBM at $\Gamma$) and Conduction Band Minimum (CBM along $\Gamma-M$) to ensure curvature convergence and reliable 2nd-derivative parabolic fitting for $m^*$.
2. **Vacuum-Level Potential Alignment (`strain/`)**:
   - In periodic boundary DFT, reference electrostatic potentials shift artificially under strain. Using `pp.x` (plot_num = 11), planar-averaged electrostatic potentials in the vacuum region ($V_{\mathrm{vac}}$) are extracted to align band edge shifts:
     $$E_{\mathrm{band}}^{\mathrm{aligned}}(\varepsilon) = E_{\mathrm{band}}(\varepsilon) - V_{\mathrm{vac}}(\varepsilon)$$
3. **Hexagonal Lattice Symmetry Fix**:
   - Strain generator scripts correctly adhere to the $120^\circ$ hexagonal basis ($v_{2x} = -\frac{a_0}{2}(1+\varepsilon)$), preventing spurious shear deformation.
4. **Harmonic Phonon Stability (`phonon/01_harmonic_phonopy/`)**:
   - Supercells ($3 \times 3 \times 1$, 27 atoms) with symmetry-reduced displacements evaluated with Quantum ESPRESSO (`tprnfor = .TRUE.`, `disk_io = 'none'`) to verify absence of imaginary phonon modes along $\Gamma-M-K-\Gamma$.

---

## Computational Environment & Dependencies

- **DFT Code**: Quantum ESPRESSO (v7.2+)
- **Transport**: BoltzTraP2 (v24.x)
- **Lattice Dynamics**: Phonopy (v4.5.0) & Phono3py (v4.2.0)
- **Python Environment**: Python 3.10+ / 3.14 with `numpy`, `scipy`, `matplotlib`, `scienceplots`, `h5py`.

---

## Contributors

- **Jauhar N.R.** (`jauhar/`)
- **Rizka Abdillah** (`rizka/`)
- **Nilam** (`nilam/`)
- **Agna** (`agna/`)
