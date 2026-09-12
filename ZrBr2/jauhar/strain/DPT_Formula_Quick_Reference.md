# Ringkasan Cepat Rumus Termoelektrik & Deformation Potential Theory (DPT)

Dokumen ini berisi **kumpulan rumus akhir (ready-to-use)** yang digunakan dalam menghitung parameter-parameter penting pada analisis termoelektrik monolayer $\mathrm{ZrXY}$ ($\mathrm{ZrFBr}$ dan $\mathrm{ZrBr}_2$).

---

## 1. Konstanta Fundamental & Satuan Standar Internasional (SI)

| Simbol | Besaran Fisika | Nilai Numerik (SI) | Satuan |
| :--- | :--- | :--- | :--- |
| $\hbar$ | Konstanta Planck Tereduksi | $1.054571817 \times 10^{-34}$ | $\text{J}\cdot\text{s}$ |
| $k_B$ | Konstanta Boltzmann | $1.380649 \times 10^{-23}$ | $\text{J/K}$ |
| $m_0$ | Massa Diam Elektron Bebas | $9.1093837 \times 10^{-31}$ | $\text{kg}$ |
| $e$ | Muatan Elementer Elektron | $1.602176634 \times 10^{-19}$ | $\text{C}$ |
| $1\ \text{Ry}$ | Konversi Satuan Rydberg QE ke Joule | $2.179872 \times 10^{-18}$ | $\text{J}$ |
| $1\ \text{\AA}$ | Konversi Angstrom ke Meter | $1.0 \times 10^{-10}$ | $\text{m}$ |

---

## 2. Parameter Struktur & Mekanik Kristal 2D

### 2.1. Luas Sel Primitif 2D ($S_0$)
Untuk kisi heksagonal 2D dengan konstanta kisi ekuilibrium $a_0$:
$$S_0 = \frac{\sqrt{3}}{2} a_0^2$$
* Satuan: $\text{m}^2$.

### 2.2. Modulus Elastisitas Dua Dimensi ($C_{2D}$)
Diperoleh dari *fitting* kuadratik energi regangan elastis $\Delta E(\varepsilon) = \frac{1}{2} K \varepsilon^2$:
$$C_{2D} = \frac{1}{S_0} \left( \frac{\partial^2 E_{\text{tot}}}{\partial \varepsilon^2} \right)_{\varepsilon=0} = \frac{K}{S_0}$$
* $\Delta E = (E_{\text{tot}}(\varepsilon) - E_0) \times 2.179872 \times 10^{-18}\ \text{J}$.
* $\varepsilon = \frac{a - a_0}{a_0}$ (regangan tanpa dimensi).
* Satuan: $\text{N/m}$ atau $\text{J/m}^2$.

---

## 3. Parameter Elektronik & Interaksi Hamburan Fonon

### 3.1. Massa Efektif Pembawa Muatan ($m^*$)
Diperoleh dari turunan kedua kelengkungan pita energi di sekitar ekstremum pita (*parabolic band fitting*):
$$\frac{1}{m^*} = \frac{1}{\hbar^2} \left( \frac{\partial^2 E(\mathbf{k})}{\partial k^2} \right)_{k=k_0}$$
* **Hole ($p$-type):** Dievaluasi pada puncak pita valensi (VBM di titik $\Gamma$) $\rightarrow m_h^*$.
* **Elektron ($n$-type):** Dievaluasi pada dasar pita konduksi (CBM) $\rightarrow m_e^*$.
* Satuan input rumus: $\text{kg}$ ($m^* = (m^*/m_0) \times m_0$).

### 3.2. Konstanta Potensial Deformasi ($E_1$)
Diperoleh dari kemiringan linier pergeseran tepi pita terhadap regangan, **wajib dikoreksi terhadap potensial elektrostatik vakum absolut ($V_{\text{vac}}$)**:
$$E_1 = \frac{\partial \left( E_{\text{edge}} - V_{\text{vac}} \right)}{\partial \varepsilon}$$
* **Hole ($p$-type):** $E_{1,h} = \frac{\partial (E_{\text{VBM}} - V_{\text{vac}})}{\partial \varepsilon}$
* **Elektron ($n$-type):** $E_{1,e} = \frac{\partial (E_{\text{CBM}} - V_{\text{vac}})}{\partial \varepsilon}$
* Satuan: $\text{eV}$. Dalam perhitungan $\tau$, konversikan ke Joule: $E_1(\text{J}) = E_1(\text{eV}) \times e$.

---

## 4. Parameter Transportasi Intrinsik DPT (Wang et al. 2023)

### 4.1. Waktu Relaksasi Hamburan Fonon Akustik ($\tau$)
Formula utama hamburan fonon akustik untuk kristal 2D (Persamaan 1 Wang et al., *J. Phys.: Condens. Matter* 35, 394001, 2023):
$$\tau(T) = \frac{2 \hbar^3 C_{2D}}{3 k_B T m^* E_1^2}$$
* Satuan SI: Detik ($\text{s}$).
* Konversi ke femtodetik: $\tau\ (\text{fs}) = \tau\ (\text{s}) \times 10^{15}$.
* Catatan: Nilai $E_1$ dimasukkan dalam satuan Joule ($E_1 \times 1.602 \times 10^{-19}\ \text{J}$), dan $m^*$ dalam satuan kg.

### 4.2. Mobilitas Pembawa Muatan ($\mu$)
Dihitung melalui relasi Drude:
$$\mu(T) = \frac{e \tau(T)}{m^*} = \frac{2 e \hbar^3 C_{2D}}{3 k_B T (m^*)^2 E_1^2}$$
* Satuan SI: $\text{m}^2/\text{Vs}$.
* Konversi praktis: $\mu\ (\text{cm}^2/\text{Vs}) = \mu\ (\text{m}^2/\text{Vs}) \times 10^4$.

---

## 5. Parameter Kinerja Termoelektrik

### 5.1. Konduktivitas Listrik Nyata ($\sigma$)
Hasil kalkulasi transport BoltzTraP2 menghasilkan koefisien tanpa waktu relaksasi $\frac{\sigma}{\tau_0}(\mu, T)$. Konduktivitas fisik dihitung dengan menyuntikkan $\tau(T)$:
$$\sigma(\mu, T) = \left[ \frac{\sigma}{\tau_0}(\mu, T) \right] \cdot \tau(T)$$
* Satuan: $\Omega^{-1}\text{m}^{-1}$ atau $\text{S/m}$.

### 5.2. Power Factor ($PF$)
Daya listrik termoelektrik per unit gradien temperatur:
$$PF(\mu, T) = S^2(\mu, T) \cdot \sigma(\mu, T)$$
* $S(\mu, T)$: Koefisien Seebeck ($\text{V/K}$).
* Satuan: $\text{W}/(\text{m}\cdot\text{K}^2)$, sering dinyatakan dalam $\text{mW}/(\text{m}\cdot\text{K}^2) = 10^{-3}\ \text{W}/(\text{m}\cdot\text{K}^2)$.

### 5.3. Konduktivitas Termal Total ($\kappa$)
Total panas yang merambat di dalam material merupakan penjumlahan dari kontribusi elektron dan fonon kisi:
$$\kappa = \kappa_e + \kappa_l$$
* **Elektronik ($\kappa_e$):** Mengikuti hukum Wiedemann-Franz:
  $$\kappa_e = L \sigma T$$
  dengan $L \approx 1.5 - 2.4 \times 10^{-8}\ \text{W}\Omega\text{K}^{-2}$ (faktor Lorentz).
* **Kisi Kristal ($\kappa_l$):** Dihitung dari interaksi fonon anharmonik, umumnya menurun seiring suhu ($\kappa_l \propto T^{-1}$).

### 5.4. Figure of Merit Termoelektrik ($ZT$)
Efisiensi konversi energi termoelektrik tanpa dimensi:
$$ZT(\mu, T) = \frac{S^2 \sigma T}{\kappa_e + \kappa_l} = \frac{PF \cdot T}{\kappa_e + \kappa_l}$$
* Semakin besar $ZT$, semakin tinggi efisiensi konversi panas menjadi listrik. Ambang batas komersial devais industri adalah $ZT \ge 1.0$.

---

## 6. Lembar Kerja Cepat (Quick Recipe)

```
[Langkah 1]  a0 (Å) ───────> S0 = (sqrt(3)/2) * a0^2
[Langkah 2]  ΔE vs ε² ─────> C2D = K / S0
[Langkah 3]  d²E/dk² ──────> m* = hbar^2 / (d²E/dk²)
[Langkah 4]  ΔE_edge vs ε ─> E1 = d(E_edge - V_vac) / dε
[Langkah 5]  C2D, m*, E1 ──> τ = (2 hbar^3 C2D) / (3 kB T m* E1^2)
[Langkah 6]  τ, m* ────────> μ = (e * τ) / m*
[Langkah 7]  BTP2 + τ ─────> σ = (σ/τ0) * τ
[Langkah 8]  S, σ ─────────> PF = S² * σ
[Langkah 9]  PF, T, κe, κl ─> ZT = (PF * T) / (κe + κl)
```
