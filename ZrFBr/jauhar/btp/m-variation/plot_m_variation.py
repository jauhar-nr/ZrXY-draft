import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    plt.style.use(['science', 'nature', 'no-latex'])
except:
    pass

m_values = [2, 16, 32]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for m in m_values:
    filename = f"interpolation_m{m}.trace"
    try:
        data = np.loadtxt(filename, skiprows=1)
    except:
        continue
    
    Ef_ry = data[:, 0]
    T = data[:, 1]
    N_uc = data[:, 2]
    S = data[:, 4] * 1e6  # uV/K
    sigma = data[:, 5]
    kappa = data[:, 7]
    
    # Filter 300K saja
    mask = (T == 300)
    
    idx_ref = np.argmin(np.abs(N_uc[mask]))
    ef_ref = Ef_ry[mask][idx_ref]
    x_eV = (Ef_ry[mask] - ef_ref) * 13.605698
    
    axes[0].plot(x_eV, S[mask], label=f'm = {m}')
    axes[1].plot(x_eV, sigma[mask], label=f'm = {m}')
    axes[2].plot(x_eV, kappa[mask], label=f'm = {m}')

axes[0].set_title('Seebeck coefficient')
axes[0].set_ylabel(r'$S$ ($\mu$V/K)')
axes[0].axhline(0, color='gray', linestyle=':', linewidth=1)

axes[1].set_title('Electrical conductivity')
axes[1].set_ylabel(r'$\sigma/\tau$ ($\Omega$ m s)$^{-1}$')

axes[2].set_title('Electronic thermal conductivity')
axes[2].set_ylabel(r'$\kappa_e/\tau$ (W m$^{-1}$ K$^{-1}$ s$^{-1}$)')

for ax in axes:
    ax.set_xlabel(r'Energy $\mu - E_f$ (eV)')
    ax.set_xlim(-1.5, 1.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
    ax.legend()

fig.suptitle('Uji Resolusi Interpolasi (-m) pada Suhu 300K', fontsize=14)
plt.tight_layout()
plt.savefig('m_variation_3panel.png', dpi=300)
