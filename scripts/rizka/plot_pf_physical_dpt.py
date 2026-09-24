#!/usr/bin/env python3
"""
plot_pf_physical_dpt.py
Physical Power Factor (PF = S^2 * sigma) and Electronic Thermal Conductivity (kappa_e)
evaluated from BoltzTraP2 transport coefficients and first-principles DPT relaxation times (tau_DPT).

Applies to 2D Janus monolayers: ZrFCl, ZrFBr, and ZrClBr.
"""

import os
import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

try:
    import scienceplots
    plt.style.use(['science', 'nature', 'no-latex'])
except Exception:
    pass

print("=== Menjalankan Perhitungan Power Factor Fisis DPT (ZrXY) ===")

# Base directory setup
script_dir = os.path.dirname(os.path.abspath(__file__))
workspace_root = os.path.abspath(os.path.join(script_dir, "../../.."))
bolz_dir = os.path.abspath(os.path.join(script_dir, "../.."))

materials = {
    'ZrFCl': {
        'title': r'Janus $\mathrm{ZrFCl}$',
        'h': 3.0909,
        'c_box': 10.0,
        'eg': 1.2744,
        'candidates': [
            os.path.join(bolz_dir, 'ZrFCl', 'interpolation.condtens'),
            os.path.join(workspace_root, 'ZrFCl', 'bolz', 'interpolation.condtens'),
        ],
        'tau': {
            300.0: {'e': 32.940e-15, 'h': 26.725e-15},
            600.0: {'e': 16.470e-15, 'h': 13.365e-15},
            900.0: {'e': 10.980e-15, 'h': 8.905e-15},
        }
    },
    'ZrFBr': {
        'title': r'Janus $\mathrm{ZrFBr}$',
        'h': 3.1985,
        'c_box': 10.0,
        'eg': 1.1938,
        'candidates': [
            os.path.join(bolz_dir, 'ZrFBr', 'interpolation.condtens'),
            os.path.join(workspace_root, 'ZrFBr', 'bolz', 'interpolation.condtens'),
        ],
        'tau': {
            300.0: {'e': 30.585e-15, 'h': 19.405e-15},
            600.0: {'e': 15.295e-15, 'h': 9.705e-15},
            900.0: {'e': 10.195e-15, 'h': 6.465e-15},
        }
    },
    'ZrClBr': {
        'title': r'Janus $\mathrm{ZrClBr}$',
        'h': 3.5432,
        'c_box': 10.0,
        'eg': 0.9030,
        'candidates': [
            os.path.join(bolz_dir, 'ZrClBr', 'interpolation.condtens'),
            os.path.join(workspace_root, 'ZrClBr', 'bolz', 'interpolation.condtens'),
        ],
        'tau': {
            300.0: {'e': 19.335e-15, 'h': 17.155e-15},
            600.0: {'e': 9.665e-15,  'h': 8.580e-15},
            900.0: {'e': 6.445e-15,  'h': 5.720e-15},
        }
    }
}

temperatures = [300.0, 600.0, 900.0]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] # 300 K (biru), 600 K (oranye), 900 K (hijau)

# Filter materials with available condtens
active_materials = {}
for mat_name, mat_info in materials.items():
    found_path = None
    for p in mat_info['candidates']:
        if os.path.exists(p) and os.path.getsize(p) > 1000:
            found_path = p
            break
    if found_path:
        mat_info['cond_file'] = found_path
        active_materials[mat_name] = mat_info
        print(f"[{mat_name}] Ditemukan file condtens: {found_path}")
    else:
        print(f"[{mat_name}] File condtens belum tersedia.")

n_mat = len(active_materials)
if n_mat == 0:
    print("Error: Tidak ada file condtens yang ditemukan!", file=sys.stderr)
    sys.exit(1)

fig, axes = plt.subplots(1, n_mat, figsize=(5.0 * n_mat, 4.5), sharey=True)
if n_mat == 1:
    axes = [axes]

optimal_summary = []

for idx, (mat_name, mat_info) in enumerate(active_materials.items()):
    ax = axes[idx]
    cond_path = mat_info['cond_file']
    print(f"Memproses {mat_name} dari {cond_path}...")
    data = np.loadtxt(cond_path)
    
    Ef_ry = data[:, 0]
    T_arr = data[:, 1]
    N_arr = data[:, 2]
    # In-plane average of xx (col 3) and yy (col 7)
    sigma_tau_arr = 0.5 * (data[:, 3] + data[:, 7])
    # In-plane average of Seebeck xx (col 12) and yy (col 16)
    S_arr = 0.5 * (data[:, 12] + data[:, 16])
    
    thick_scale = mat_info['c_box'] / mat_info['h']
    eg = mat_info['eg']
    v_min = -eg / 2.0
    c_max = eg / 2.0
    
    # Shade band gap
    ax.axvspan(v_min, c_max, color='#f0f0f0', alpha=0.85, zorder=0)
    ax.axvline(0, color='gray', linestyle='--', lw=0.6, alpha=0.7)
    ax.text(0.0, 28, f'$E_g \\approx {eg:.2f}\\ \\mathrm{{eV}}$', 
            ha='center', va='center', fontsize=8, color='dimgray',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.85, edgecolor='lightgray'))
    
    for t_idx, temp in enumerate(temperatures):
        mask_T = T_arr == temp
        if not np.any(mask_T):
            print(f"[{mat_name}] Suhu {temp} K tidak ditemukan dalam data, dilewati.")
            continue
        ef_T = Ef_ry[mask_T]
        n_T = N_arr[mask_T]
        sigma_tau = sigma_tau_arr[mask_T] * thick_scale
        S = S_arr[mask_T]
        
        # Cari tingkat Fermi intrinsik (N terdekat dengan 0 / netral)
        idx_int = np.argmin(np.abs(n_T))
        ef_ref = ef_T[idx_int]
        
        x_eV = (ef_T - ef_ref) * 13.605698
        sort_idx = np.argsort(x_eV)
        x_sort = x_eV[sort_idx]
        
        range_mask = (x_sort >= -1.5) & (x_sort <= 1.5)
        x_plot = x_sort[range_mask]
        
        sig_plot = sigma_tau[sort_idx][range_mask]
        s_plot = S[sort_idx][range_mask]
        
        tau_e = mat_info['tau'][temp]['e']
        tau_h = mat_info['tau'][temp]['h']
        
        # Terapkan tau_h untuk mu - Ef < 0 (hole), dan tau_e untuk mu - Ef > 0 (elektron)
        tau_array = np.where(x_plot < 0, tau_h, tau_e)
        pf_phys = (s_plot**2) * sig_plot * tau_array * 1e3 # mW / (m K^2)
        
        c = colors[t_idx]
        ax.plot(x_plot, pf_phys, color=c, lw=1.6, label=f'{int(temp)} K')
        
        # Puncak p-type (mu - Ef < v_min) dan n-type (mu - Ef > c_max)
        mask_p = (x_plot < -0.1)
        mask_n = (x_plot > 0.1)
        
        if np.any(mask_p):
            max_idx_p = np.argmax(pf_phys[mask_p])
            peak_p = pf_phys[mask_p][max_idx_p]
            mu_opt_p = x_plot[mask_p][max_idx_p]
            ax.scatter([mu_opt_p], [peak_p], color=c, edgecolors='black', s=35, marker='s', zorder=5)
        else:
            peak_p, mu_opt_p = 0.0, 0.0
            
        if np.any(mask_n):
            max_idx_n = np.argmax(pf_phys[mask_n])
            peak_n = pf_phys[mask_n][max_idx_n]
            mu_opt_n = x_plot[mask_n][max_idx_n]
            ax.scatter([mu_opt_n], [peak_n], color=c, edgecolors='black', s=35, marker='o', zorder=5)
        else:
            peak_n, mu_opt_n = 0.0, 0.0
            
        optimal_summary.append({
            'mat': mat_name,
            'T': int(temp),
            'PF_max_p': peak_p,
            'mu_opt_p': mu_opt_p,
            'PF_max_n': peak_n,
            'mu_opt_n': mu_opt_n
        })
        
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(0, 35)
    ax.grid(True, linestyle=':', lw=0.4, alpha=0.6)
    ax.set_xlabel(r'$\mu - E_F$ (eV)', fontsize=10)
    ax.set_title(f'({chr(97+idx)}) {mat_info["title"]}', fontsize=11, fontweight='bold')
    
    if idx == 0:
        ax.set_ylabel(r'Physical Power Factor $PF$ ($\mathrm{mW}\,\mathrm{m}^{-1}\mathrm{K}^{-2}$)', fontsize=10)
        ax.legend(frameon=True, fontsize=8, loc='upper right', framealpha=0.9)
    ax.text(-0.95, 23, r'$\leftarrow p$-type (Hole)', fontsize=8, color='#444444', ha='center')
    ax.text(1.05, 23, r'$n$-type (Elec) $\rightarrow$', fontsize=8, color='#444444', ha='center')

plt.tight_layout()
out_dir = os.path.join(script_dir, "output")
os.makedirs(out_dir, exist_ok=True)
fig_path1 = os.path.join(script_dir, "Fig_PF_Physical_DPT.png")
fig_path2 = os.path.join(workspace_root, "artifacts", "Fig_PF_Physical_DPT.png")

plt.savefig(fig_path1, dpi=400, bbox_inches='tight')
plt.savefig(fig_path2, dpi=400, bbox_inches='tight')
print(f"[SUKSES] Gambar disimpan di:\n  - {fig_path1}\n  - {fig_path2}")

# Print and save summary table
table_md_path = os.path.join(script_dir, "Table_PF_Benchmarking.md")
table_csv_path = os.path.join(script_dir, "Table_PF_Benchmarking.csv")
art_table_path = os.path.join(workspace_root, "artifacts", "Table_PF_Benchmarking.md")

lines = [
    "# Tabel Ringkasan Power Factor Fisis DPT Monolayer Janus ZrXY",
    "",
    "| Material | T (K) | PF_max (p-type) [mW/mK^2] | (mu-E_F)_opt (p) [eV] | PF_max (n-type) [mW/mK^2] | (mu-E_F)_opt (n) [eV] |",
    "| :--- | :---: | :---: | :---: | :---: | :---: |",
]

csv_lines = ["Material,T_K,PF_max_p_mW_mK2,mu_opt_p_eV,PF_max_n_mW_mK2,mu_opt_n_eV\n"]

for r in optimal_summary:
    lines.append(f"| **{r['mat']}** | {r['T']} | **{r['PF_max_p']:.2f}** | {r['mu_opt_p']:.3f} | **{r['PF_max_n']:.2f}** | {r['mu_opt_n']:.3f} |")
    csv_lines.append(f"{r['mat']},{r['T']},{r['PF_max_p']:.4f},{r['mu_opt_p']:.4f},{r['PF_max_n']:.4f},{r['mu_opt_n']:.4f}\n")

with open(table_md_path, 'w') as f:
    f.write("\n".join(lines) + "\n")
with open(art_table_path, 'w') as f:
    f.write("\n".join(lines) + "\n")
with open(table_csv_path, 'w') as f:
    f.writelines(csv_lines)

print(f"[SUKSES] Tabel ringkasan disimpan di:\n  - {table_md_path}\n  - {art_table_path}\n  - {table_csv_path}")

print("\n" + "="*85)
print("RINGKASAN PUNCAK POWER FACTOR FISIS DPT:")
print("="*85)
for r in optimal_summary:
    print(f"{r['mat']:<8} | {r['T']:<5} K | PF_max(p): {r['PF_max_p']:<6.2f} at {r['mu_opt_p']:<6.3f} eV | PF_max(n): {r['PF_max_n']:<6.2f} at {r['mu_opt_n']:<6.3f} eV")
print("="*85)
