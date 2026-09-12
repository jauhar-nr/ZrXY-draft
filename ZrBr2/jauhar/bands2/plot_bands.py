import numpy as np
import matplotlib.pyplot as plt
import os
import re

try:
    import scienceplots
    plt.style.use(['science', 'nature', 'no-latex'])
except Exception:
    pass

print("Membaca data Super-Dense Band Structure untuk ZrBr2...")

# 1. Baca Fermi Energy dari scf.out
ef = 0.0
try:
    with open('scf.out', 'r') as f:
        content = f.read()
    m1 = re.search(r'Fermi energy is\s+([+-]?\d+\.\d+)', content)
    m2 = re.search(r'highest occupied.*:\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)', content)
    if m1:
        ef = float(m1.group(1))
    elif m2:
        ef = (float(m2.group(1)) + float(m2.group(2))) / 2.0
except Exception as err:
    print(f"[WARN] Gagal membaca Fermi level dari scf.out: {err}")

# 2. Baca High-Symmetry Ticks dari bands.out
ticks = []
labels = [r'$\Gamma$', 'M', 'K', r'$\Gamma$']
try:
    with open('bands.out', 'r') as f:
        for line in f:
            if "high-symmetry point:" in line:
                x_coord = float(line.split()[-1])
                ticks.append(x_coord)
except Exception:
    ticks = [0.0, 0.5774, 0.9107, 1.5774]

if not ticks:
    ticks = [0.0, 0.5774, 0.9107, 1.5774]

# 3. Baca Data Pita Energi dari bands.dat.gnu
bands = []
current_band = []
with open('bands.dat.gnu', 'r') as f:
    for line in f:
        if line.strip() == "":
            if len(current_band) > 0:
                bands.append(np.array(current_band))
                current_band = []
        else:
            parts = line.split()
            current_band.append([float(parts[0]), float(parts[1])])
if len(current_band) > 0:
    bands.append(np.array(current_band))

print(f"Total pita terdeteksi: {len(bands)}, titik k per pita: {len(bands[0])}")

# 4. Analisis Titik Ekstrem (VBM & CBM) dan Massa Efektif
# Jumlah elektron = 26.00 -> occupied bands: 1 sampai 13 (index 0 sampai 12)
# VBM: band 13 (index 12); CBM: band 14 (index 13)
vb = bands[12]
cb = bands[13]

vbm_idx = np.argmax(vb[:, 1])
cbm_idx = np.argmin(cb[:, 1])

vbm_E = vb[vbm_idx, 1]
cbm_E = cb[cbm_idx, 1]
band_gap = cbm_E - vbm_E

vbm_x = vb[vbm_idx, 0]
cbm_x = cb[cbm_idx, 0]

# Konversi koordinat k ke satuan Angstrom^-1
a_lat = 3.560487104 # Angstrom
scale_k = 2.0 * np.pi / a_lat # 1/A
HBAR2_OVER_M0 = 7.619964 # eV * A^2

# Fitting Parabolik di sekitar VBM (+/- 4 titik)
k_vbm_fit = (vb[vbm_idx-4:vbm_idx+5, 0] - vbm_x) * scale_k
E_vbm_fit = vb[vbm_idx-4:vbm_idx+5, 1]
p_vb = np.polyfit(k_vbm_fit, E_vbm_fit, 2)
m_h = -HBAR2_OVER_M0 / (2.0 * p_vb[0])

# Fitting Parabolik di sekitar CBM (+/- 4 titik)
k_cbm_fit = (cb[cbm_idx-4:cbm_idx+5, 0] - cbm_x) * scale_k
E_cbm_fit = cb[cbm_idx-4:cbm_idx+5, 1]
p_cb = np.polyfit(k_cbm_fit, E_cbm_fit, 2)
m_e = HBAR2_OVER_M0 / (2.0 * p_cb[0])

summary_text = f"""==================================================
Hasil Analisis Super-Dense Band Structure: ZrBr2
==================================================
Kisi a               : {a_lat:.9f} Angstrom
Fermi Energy         : {ef:.4f} eV
VBM Energy           : {vbm_E:.4f} eV (pada x = {vbm_x:.4f})
CBM Energy           : {cbm_E:.4f} eV (pada x = {cbm_x:.4f})
Band Gap (Indirect)  : {band_gap:.4f} eV
Massa Efektif Hole   : {m_h:.4f} m_0
Massa Efektif Elec   : {m_e:.4f} m_0
==================================================
"""
print(summary_text)
with open('effective_mass.txt', 'w') as f_out:
    f_out.write(summary_text)

# 5. Plotting Band Structure Resolusi Tinggi
fig, ax = plt.subplots(figsize=(4.2, 5.2))

for b in bands:
    ax.plot(b[:, 0], b[:, 1] - ef, color='black', lw=1.1)

# Highlight VBM dan CBM
ax.scatter([vbm_x], [vbm_E - ef], color='blue', s=25, zorder=5, label=f'VBM ($m_h^*={m_h:.2f}m_0$)')
ax.scatter([cbm_x], [cbm_E - ef], color='red', s=25, zorder=5, label=f'CBM ($m_e^*={m_e:.2f}m_0$)')

for x in ticks:
    ax.axvline(x=x, color='gray', lw=0.5, ls='--')

ax.axhline(y=0, color='red', lw=0.8, ls='--', alpha=0.7)

ax.set_xlim(ticks[0], ticks[-1])
ax.set_ylim(-5, 5)

ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel(r'Energy - $E_F$ (eV)', fontsize=11)
ax.set_title(r'$\mathrm{ZrBr}_2$ Monolayer ($E_g = ' + f'{band_gap:.2f}' + r'\ \mathrm{eV}$)', fontsize=11)
ax.legend(frameon=True, fontsize=8, loc='upper right')

plt.tight_layout()
plt.savefig('BandStructure_SuperDense.png', dpi=600, bbox_inches='tight')
print("Gambar disimpan sebagai BandStructure_SuperDense.png")
