import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

try:
    import scienceplots
    plt.style.use(['science', 'nature', 'no-latex'])
except Exception:
    pass

# Konstanta Fundamental
HBAR = 1.054571817e-34       # J s
KB = 1.380649e-23           # J/K
M0 = 9.1093837e-31          # kg
E_CHG = 1.602176634e-19     # C
RY_TO_J = 2.179872e-18      # J/Ry
ANG_TO_M = 1.0e-10          # m

a0 = 3.343942779 # Angstrom
S0_m2 = (np.sqrt(3.0) / 2.0) * ((a0 * ANG_TO_M)**2) # m^2

# 1. Baca massa efektif dari bands2 jika tersedia
m_h_default = 0.729
m_e_default = 1.276
eff_file = "../bands2/effective_mass.txt"
if os.path.exists(eff_file):
    try:
        with open(eff_file, 'r') as f:
            content = f.read()
        mh_match = re.search(r'Massa Efektif Hole\s+:\s+([0-9\.]+)', content)
        me_match = re.search(r'Massa Efektif Elec\s+:\s+([0-9\.]+)', content)
        if mh_match: m_h_default = float(mh_match.group(1))
        if me_match: m_e_default = float(me_match.group(1))
    except Exception:
        pass

def parse_relax_out(filepath):
    """Membaca energi total akhir dan tingkat energi VBM/CBM"""
    e_tot = None
    vbm = None
    cbm = None
    if not os.path.exists(filepath):
        return None, None, None
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if "!" in line and "total energy" in line:
            parts = line.split()
            e_tot = float(parts[-2]) # Ry
        elif "highest occupied, lowest unoccupied level" in line:
            parts = line.split()
            vbm = float(parts[-2])
            cbm = float(parts[-1])
        elif "highest occupied level" in line and vbm is None:
            parts = line.split()
            vbm = float(parts[-1])
    return e_tot, vbm, cbm

def parse_vacuum_level(potpath):
    """Membaca potensial vakum dari pot.dat"""
    if not os.path.exists(potpath):
        return 0.0
    try:
        data = []
        with open(potpath, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        data.append([float(parts[0]), float(parts[1])])
                    except ValueError:
                        pass
        if len(data) > 0:
            arr = np.array(data)
            # Ambil daerah vakum (dekat z = 0 atau z = c/2)
            # Slab terletak di sekitar tengah (z ~ 10 A), vakum di z ~ 0 - 4 A
            vac_mask = (arr[:, 0] < 4.0) | (arr[:, 0] > 16.0)
            if np.any(vac_mask):
                return np.mean(arr[vac_mask, 1])
            return arr[0, 1]
    except Exception:
        pass
    return 0.0

directions = ['dir_x', 'dir_y']
results = {}

fig, axes = plt.subplots(2, 2, figsize=(9, 7))

for d_idx, d in enumerate(directions):
    subs = sorted(glob.glob(f"{d}/eps_*"))
    eps_list = []
    e_tot_list = []
    vbm_list = []
    cbm_list = []
    
    for sub in subs:
        m = re.search(r'eps_([+-][0-9\.]+)', sub)
        if not m: continue
        eps = float(m.group(1))
        
        relax_out = os.path.join(sub, "relax.out")
        pot_dat = os.path.join(sub, "pot.dat")
        
        etot, vbm, cbm = parse_relax_out(relax_out)
        v_vac = parse_vacuum_level(pot_dat)
        
        if etot is not None and vbm is not None and cbm is not None:
            eps_list.append(eps)
            e_tot_list.append(etot)
            # Level relatif terhadap vakum
            vbm_list.append(vbm - v_vac)
            cbm_list.append(cbm - v_vac)
            
    if len(eps_list) < 3:
        print(f"[WARN] Belum cukup data kalkulasi di {d} (ditemukan {len(eps_list)} titik)")
        continue
        
    eps_arr = np.array(eps_list)
    etot_arr = np.array(e_tot_list)
    vbm_arr = np.array(vbm_list)
    cbm_arr = np.array(cbm_list)
    
    # 1. Fit Modulus Elastis C_2D: Delta E = 1/2 * K * eps^2 (Joule)
    e_rel_J = (etot_arr - np.min(etot_arr)) * RY_TO_J
    p_elast = np.polyfit(eps_arr, e_rel_J, 2)
    # p_elast[0] = 1/2 * K -> K = 2 * p_elast[0]
    # C_2D = K / S0 (J/m^2 = N/m)
    C_2D = (2.0 * p_elast[0]) / S0_m2
    
    # 2. Fit Potensial Deformasi E1: Delta E_edge = E1 * eps (eV)
    p_vbm = np.polyfit(eps_arr, vbm_arr, 1)
    p_cbm = np.polyfit(eps_arr, cbm_arr, 1)
    E1_h = p_vbm[0] # eV
    E1_e = p_cbm[0] # eV
    
    # 3. Hitung Waktu Relaksasi tau dan Mobilitas mu (Persamaan 1 Zr2Cl4)
    # tau = (2 * hbar^3 * C) / (3 * kB * T * m* * E1^2)
    temps = [300.0, 600.0, 900.0]
    
    def calc_tau_mu(m_eff_rel, E1_val, T_val):
        m_eff = m_eff_rel * M0
        E1_J = E1_val * E_CHG
        tau_s = (2.0 * (HBAR**3) * C_2D) / (3.0 * KB * T_val * m_eff * (E1_J**2))
        mu_m2Vs = (E_CHG * tau_s) / m_eff
        mu_cm2Vs = mu_m2Vs * 1e4
        return tau_s * 1e15, mu_cm2Vs # fs, cm^2/Vs
        
    results[d] = {
        'C_2D': C_2D,
        'E1_h': E1_h,
        'E1_e': E1_e,
        'm_h': m_h_default,
        'm_e': m_e_default,
        'temps': {}
    }
    
    for T in temps:
        tau_h, mu_h = calc_tau_mu(m_h_default, E1_h, T)
        tau_e, mu_e = calc_tau_mu(m_e_default, E1_e, T)
        results[d]['temps'][int(T)] = {
            'tau_h': tau_h, 'mu_h': mu_h,
            'tau_e': tau_e, 'mu_e': mu_e
        }
        
    # Plotting
    # Baris 1: Energi vs Strain
    ax_e = axes[0, d_idx]
    eps_fit = np.linspace(np.min(eps_arr), np.max(eps_arr), 100)
    ax_e.scatter(eps_arr * 100, e_rel_J * 1e18, color='black', s=25, label='DFT')
    ax_e.plot(eps_fit * 100, np.polyval(p_elast, eps_fit) * 1e18, color='crimson', lw=1.2,
              label=f'Parabolic Fit\n$C_{{2D}} = {C_2D:.2f}$ N/m')
    ax_e.set_xlabel(r'Strain $\varepsilon$ (%)')
    ax_e.set_ylabel(r'$\Delta E_{\mathrm{tot}}$ ($10^{-18}$ J)')
    ax_e.set_title(f'ZrFBr Elastic Energy ({d})')
    ax_e.legend(frameon=True, fontsize=8)
    
    # Baris 2: Band Edges vs Strain
    ax_b = axes[1, d_idx]
    ax_b.scatter(eps_arr * 100, vbm_arr, color='blue', s=25, label='VBM')
    ax_b.plot(eps_fit * 100, np.polyval(p_vbm, eps_fit), color='blue', ls='--', lw=1.1,
              label=f'VBM Fit ($E_{{1,h}}={E1_h:.2f}$ eV)')
    ax_b.scatter(eps_arr * 100, cbm_arr, color='red', s=25, label='CBM')
    ax_b.plot(eps_fit * 100, np.polyval(p_cbm, eps_fit), color='red', ls='--', lw=1.1,
              label=f'CBM Fit ($E_{{1,e}}={E1_e:.2f}$ eV)')
    ax_b.set_xlabel(r'Strain $\varepsilon$ (%)')
    ax_b.set_ylabel('Band Edge vs Vacuum (eV)')
    ax_b.set_title(f'ZrFBr Deformation Potential ({d})')
    ax_b.legend(frameon=True, fontsize=8)

plt.tight_layout()
plt.savefig('Strain_Fitting.png', dpi=600, bbox_inches='tight')
print("Plot disimpan sebagai Strain_Fitting.png")

# Tulis Hasil Lengkap ke File Teks
with open('DPT_results.txt', 'w') as f_out:
    f_out.write("======================================================================\n")
    f_out.write("HASIL ANALISIS DEFORMATION POTENTIAL THEORY (DPT): ZrFBr\n")
    f_out.write("Formula: Persamaan (1) Wang et al. JPCM 35 (2023) 394001\n")
    f_out.write("======================================================================\n\n")
    
    for d, res in results.items():
        f_out.write(f"--- Arah: {d} ---\n")
        f_out.write(f"Modulus Elastis C_2D : {res['C_2D']:.2f} N/m (J/m^2)\n")
        f_out.write(f"Massa Efektif Hole   : {res['m_h']:.4f} m_0\n")
        f_out.write(f"Massa Efektif Elec   : {res['m_e']:.4f} m_0\n")
        f_out.write(f"DP Constant Hole E1  : {res['E1_h']:.3f} eV\n")
        f_out.write(f"DP Constant Elec E1  : {res['E1_e']:.3f} eV\n\n")
        f_out.write("Suhu (K) | Tipe Carrier | Waktu Relaksasi tau (fs) | Mobilitas mu (cm^2/Vs)\n")
        f_out.write("----------------------------------------------------------------------\n")
        for T, vals in res['temps'].items():
            f_out.write(f"{T:8d} | Electron (n) | {vals['tau_e']:24.2f} | {vals['mu_e']:20.2f}\n")
            f_out.write(f"{T:8d} | Hole (p)     | {vals['tau_h']:24.2f} | {vals['mu_h']:20.2f}\n")
        f_out.write("\n")

print("Analisis selesai! Rangkuman disimpan di DPT_results.txt")
