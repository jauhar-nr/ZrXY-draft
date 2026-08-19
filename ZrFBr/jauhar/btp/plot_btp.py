import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science', 'nature', 'no-latex'])

print("Membaca data dari interpolation.trace...")
data = np.loadtxt('interpolation.trace')

Ef_ry = data[:, 0]
T = data[:, 1]
N_uc = data[:, 2]
S = data[:, 4] * 1e6  # uV/K
sigma = data[:, 5] / 1e18 # 10^18 /(ohm m s)
kappa_e = data[:, 7] / 1e14 # 10^14 W/(m K s)

temperatures = np.unique(T)

fig, axs = plt.subplots(1, 3, figsize=(11, 3))

for temp in temperatures:
    mask = T == temp
    ef = Ef_ry[mask]
    n = N_uc[mask]
    s = S[mask]
    sig = sigma[mask]
    kap = kappa_e[mask]
    
    # Intrinsic Fermi level = energy di mana doping (N) mendekati 0 untuk suhu ini
    idx_int = np.argmin(np.abs(n))
    ef_int = ef[idx_int]
    
    # 1 Ry = 13.605698 eV
    x_eV = (ef - ef_int) * 13.605698
    
    axs[0].plot(x_eV, s, label=f'{int(temp)} K', lw=1.2)
    axs[1].plot(x_eV, sig, lw=1.2)
    axs[2].plot(x_eV, kap, lw=1.2)

# Styling Plot 1: Seebeck
axs[0].set_xlim(-1.5, 1.5)
axs[0].set_xlabel(r'$\mu - E_F$ (eV)')
axs[0].set_ylabel(r'Seebeck Coefficient ($\mu$V/K)')
axs[0].axvline(0, color='black', linestyle='--', lw=0.5, alpha=0.5)
axs[0].axhline(0, color='black', linestyle='--', lw=0.5, alpha=0.5)
axs[0].legend(frameon=False, fontsize=8, loc='best')

# Styling Plot 2: Electrical Conductivity
axs[1].set_xlim(-1.5, 1.5)
# axs[1].set_ylim(bottom=0)
axs[1].set_ylim(0, 250)
axs[1].set_xlabel(r'$\mu - E_F$ (eV)')
axs[1].set_ylabel(r'$\sigma / \tau$ ($10^{18}$ $\Omega^{-1}$m$^{-1}$s$^{-1}$)')
axs[1].axvline(0, color='black', linestyle='--', lw=0.5, alpha=0.5)

# Styling Plot 3: Electronic Thermal Conductivity (kappa_e)
axs[2].set_xlim(-1.5, 1.5)
# axs[2].set_ylim(bottom=0)
axs[2].set_ylim(0, 45)
axs[2].set_xlabel(r'$\mu - E_F$ (eV)')
axs[2].set_ylabel(r'$\kappa_e / \tau$ ($10^{14}$ W m$^{-1}$ K$^{-1}$ s$^{-1}$)')
axs[2].axvline(0, color='black', linestyle='--', lw=0.5, alpha=0.5)

plt.tight_layout()
plt.savefig('Thermoelectric_Properties.png', dpi=400, bbox_inches='tight')
print("Selesai! Grafik disimpan sebagai Thermoelectric_Properties.png")
