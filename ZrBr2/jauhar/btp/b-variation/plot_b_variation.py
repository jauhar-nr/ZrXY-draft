import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    plt.style.use(['science', 'nature', 'no-latex'])
except:
    pass

b_values = [1000, 5000, 10000, 50000]

fig, ax = plt.subplots(figsize=(8, 6))

for b in b_values:
    filename = f"interpolation_b{b}.trace"
    try:
        data = np.loadtxt(filename, skiprows=1)
    except:
        print(f"File {filename} tidak ditemukan, melompati...")
        continue
    
    Ef_ry = data[:, 0]
    T = data[:, 1]
    N_uc = data[:, 2]
    S = data[:, 4] * 1e6  # uV/K
    
    # Ambil data suhu 300K
    mask_300 = (T == 300)
    
    # Cari Refined Fermi Level (N = 0)
    idx_ref = np.argmin(np.abs(N_uc[mask_300]))
    ef_ref = Ef_ry[mask_300][idx_ref]
    
    # Konversi sumbu X ke eV
    x_eV = (Ef_ry[mask_300] - ef_ref) * 13.605698
    
    ax.plot(x_eV, S[mask_300], label=f'b = {b}')

ax.set_xlabel(r'Energy $(E - E_f)$ (eV)')
ax.set_ylabel(r'Seebeck Coefficient ($\mu$V/K)')
ax.set_xlim(-2, 2)
ax.legend()
ax.set_title('Uji Variasi Resolusi Bins (b) pada 300K')

plt.tight_layout()
plt.savefig('b_variation_plot.png', dpi=300)
print("Plot berhasil disimpan sebagai b_variation_plot.png")
