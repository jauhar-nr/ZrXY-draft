# Teori Potensial Deformasi (Deformation Potential Theory - DPT): Penurunan Matematis, Metodologi Komputasi DFT, dan Tabel Hasil Monolayer $\mathrm{ZrFBr}$ & $\mathrm{ZrBr}_2$

Dokumen ini menyajikan secara lengkap dan mandiri:
1. **Penurunan Matematis Lengkap** formula waktu relaksasi ($\tau$) dan mobilitas pembawa muatan ($\mu$) berbasis *Deformation Potential Theory* (DPT) untuk material dua dimensi (2D).
2. **Metodologi Komputasi DFT** untuk mengekstraksi seluruh parameter fisis ($C_{2D}$, $E_1$, $m^*$, serta koreksi potensial elektrostatik vakum absolut).
3. **Tabel Hasil Kuantitatif Lengkap** untuk monolayer Janus $\mathrm{ZrFBr}$ dan monolayer parent $\mathrm{ZrBr}_2$ pada variasi temperatur $300\ \text{K}$, $600\ \text{K}$, dan $900\ \text{K}$.

---

## 1. Pendahuluan dan Konsep Dasar DPT

Dalam kristal semikonduktor murni tanpa cacat struktural (*defect-free*), mekanisme hamburan intrinsik utama yang membatasi transportasi pembawa muatan pada suhu ruang ($T \ge 300\ \text{K}$) adalah **hamburan oleh fonon akustik** (*acoustic phonon scattering*). 

Teori Potensial Deformasi (*Deformation Potential Theory* - DPT) yang dirintis oleh **John Bardeen dan William Shockley (1950)** didasarkan pada prinsip fisis bahwa gelombang akustik berpanjang gelombang panjang ($q \to 0$) menyebabkan modulasi regangan volume lokal kisi secara berkala. Modulasi kisi ini menggeser tingkat energi pita konduksi (CBM) dan pita valensi (VBM), sehingga bertindak sebagai potensial perturbasi hamburan bagi elektron dan hole.

Pada material monolayer 2D, teori ini diperluas (Lang et al., *Phys. Rev. B* 94, 235306, 2016; Wang et al., *J. Phys.: Condens. Matter* 35, 394001, 2023) dengan mempertimbangkan:
- Modulus elastisitas kisi dua dimensi $C_{2D}$ (satuan $\text{N/m}$ atau $\text{J/m}^2$),
- Konstanta potensial deformasi $E_1$ (satuan $\text{eV}$),
- Massa efektif pembawa muatan $m^*$ (satuan $m_0$),
- Kerapatan keadaan (*density of states* - DOS) kuantum 2D yang konstan terhadap energi.

---

## 2. Penurunan Matematis Rumus DPT (Step-by-Step Derivation)

### 2.1. Hamiltonian Interaksi Elektron-Fonon Akustik

Gelombang fonon akustik longitudinal (*Longitudinal Acoustic* - LA) merambat melalui kisi kristal dengan menghasilkan perpindahan atomik lokal $\mathbf{u}(\mathbf{r})$. Perubahan volume relatif (dilatasi fraksional atau regangan lokal) didefinisikan oleh divergensi medan perpindahan:
$$\Delta(\mathbf{r}) = \nabla \cdot \mathbf{u}(\mathbf{r})$$

Pergeseran energi tepi pita akibat dilatasi lokal ini dimodelkan melalui konstanta potensial deformasi $E_1$:
$$H_{\text{ep}}(\mathbf{r}) = E_1 \Delta(\mathbf{r}) = E_1 \nabla \cdot \mathbf{u}(\mathbf{r})$$
di mana:
- $E_1 = \frac{\partial E_{\text{edge}}}{\partial \varepsilon}$ adalah konstanta potensial deformasi ($\text{eV}$).

---

### 2.2. Kuantisasi Medan Perpindahan Kisi 2D

Dalam formalisme mekanika kuantum, operator perpindahan atomik 2D $\mathbf{u}(\mathbf{r})$ dikuantisasi ke dalam penjumlahan modus fonon dengan vektor gelombang $\mathbf{q}$ dan polarisasi $\lambda$:
$$\mathbf{u}(\mathbf{r}) = \sum_{\mathbf{q}} \sqrt{\frac{\hbar}{2 \rho_{2D} A \omega_{\mathbf{q}}}} \mathbf{e}_{\mathbf{q}} \left( a_{\mathbf{q}} e^{i \mathbf{q} \cdot \mathbf{r}} + a_{\mathbf{q}}^\dagger e^{-i \mathbf{q} \cdot \mathbf{r}} \right)$$

di mana:
- $\rho_{2D} = \frac{M_{\text{cell}}}{S_0}$ adalah densitas massa kristal 2D ($\text{kg/m}^2$),
- $A$ adalah luas kristal makroskopik untuk normalisasi,
- $\omega_{\mathbf{q}} = v_s q$ adalah frekuensi fonon akustik dengan kecepatan suara $v_s$,
- $\mathbf{e}_{\mathbf{q}}$ adalah vektor polarisasi satuan mode fonon,
- $a_{\mathbf{q}}$ dan $a_{\mathbf{q}}^\dagger$ adalah operator anihilasi dan kreasi fonon.

Untuk fonon akustik longitudinal (LA), vektor perpindahan sejajar dengan arah propagasi gelombang, sehingga $\mathbf{e}_{\mathbf{q}} = \frac{\mathbf{q}}{q}$. Divergensi dari $\mathbf{u}(\mathbf{r})$ menjadi:
$$\nabla \cdot \mathbf{u}(\mathbf{r}) = \sum_{\mathbf{q}} \sqrt{\frac{\hbar}{2 \rho_{2D} A \omega_{\mathbf{q}}}} (i \mathbf{q} \cdot \mathbf{e}_{\mathbf{q}}) \left( a_{\mathbf{q}} e^{i \mathbf{q} \cdot \mathbf{r}} - a_{\mathbf{q}}^\dagger e^{-i \mathbf{q} \cdot \mathbf{r}} \right)$$

Karena $\mathbf{q} \cdot \mathbf{e}_{\mathbf{q}} = q$, maka:
$$H_{\text{ep}}(\mathbf{r}) = i E_1 \sum_{\mathbf{q}} \sqrt{\frac{\hbar q^2}{2 \rho_{2D} A v_s q}} \left( a_{\mathbf{q}} e^{i \mathbf{q} \cdot \mathbf{r}} - a_{\mathbf{q}}^\dagger e^{-i \mathbf{q} \cdot \mathbf{r}} \right)$$

---

### 2.3. Elemen Matriks Hamburan dan Batas Temperatur Tinggi

Kita tinjau proses transisi pembawa muatan dari keadaan awal Bloch $|\mathbf{k}\rangle = \frac{1}{\sqrt{A}} e^{i \mathbf{k} \cdot \mathbf{r}}$ ke keadaan akhir $|\mathbf{k}'\rangle = \frac{1}{\sqrt{A}} e^{i \mathbf{k}' \cdot \mathbf{r}}$ akibat interaksi dengan fonon $\mathbf{q}$.

Berdasarkan aturan integrasi gelombang bidang:
$$\frac{1}{A} \int e^{i (\mathbf{k} \pm \mathbf{q} - \mathbf{k}') \cdot \mathbf{r}} d^2\mathbf{r} = \delta_{\mathbf{k}', \mathbf{k} \pm \mathbf{q}}$$
Momentum kristal terkonservasi: $\mathbf{k}' = \mathbf{k} + \mathbf{q}$ (absorpsi) atau $\mathbf{k}' = \mathbf{k} - \mathbf{q}$ (emisi).

Kuadrat elemen matriks untuk kedua proses:
- **Absorpsi Fonon:**
  $$|M_{\text{abs}}|^2 = |\langle \mathbf{k} + \mathbf{q}, N_{\mathbf{q}} - 1 | H_{\text{ep}} | \mathbf{k}, N_{\mathbf{q}} \rangle|^2 = \frac{\hbar q E_1^2}{2 \rho_{2D} A v_s} N_{\mathbf{q}}$$
- **Emisi Fonon:**
  $$|M_{\text{emi}}|^2 = |\langle \mathbf{k} - \mathbf{q}, N_{\mathbf{q}} + 1 | H_{\text{ep}} | \mathbf{k}, N_{\mathbf{q}} \rangle|^2 = \frac{\hbar q E_1^2}{2 \rho_{2D} A v_s} (N_{\mathbf{q}} + 1)$$

Di mana $N_{\mathbf{q}}$ adalah fungsi distribusi Bose-Einstein:
$$N_{\mathbf{q}} = \frac{1}{\exp\left(\frac{\hbar \omega_{\mathbf{q}}}{k_B T}\right) - 1}$$

Pada temperatur ruang dan temperatur tinggi ($T \ge 300\ \text{K}$), energi fonon akustik berpanjang gelombang panjang jauh lebih kecil daripada energi termal ($\hbar \omega_{\mathbf{q}} = \hbar v_s q \ll k_B T$). Dengan ekspansi deret Taylor:
$$N_{\mathbf{q}} \approx \frac{k_B T}{\hbar \omega_{\mathbf{q}}} = \frac{k_B T}{\hbar v_s q} \gg 1$$

Sehingga probabilitas absorpsi dan emisi menjadi setara ($N_{\mathbf{q}} + 1 \approx N_{\mathbf{q}}$), dan total kuadrat elemen matriks adalah:
$$|M_{\text{tot}}|^2 = |M_{\text{abs}}|^2 + |M_{\text{emi}}|^2 = 2 \times \frac{\hbar q E_1^2}{2 \rho_{2D} A v_s} \left(\frac{k_B T}{\hbar v_s q}\right) = \frac{k_B T E_1^2}{\rho_{2D} A v_s^2}$$

Dalam mekanika elastisitas 2D, modulus elastisitas in-plane $C_{2D}$ berhubungan langsung dengan kerapatan massa 2D dan kecepatan suara fonon longitudinal:
$$C_{2D} = \rho_{2D} v_s^2$$

Substitusi hubungan ini menghasilkan bentuk matriks hamburan yang luar biasa elegan:
$$|M_{\text{tot}}|^2 = \frac{k_B T E_1^2}{A C_{2D}}$$

> **Catatan Fisis Penting:** Kuadrat elemen matriks hamburan $|M_{\text{tot}}|^2$ bersifat **independen terhadap besar maupun sudut momentum transfer $\mathbf{q}$**. Ini berarti hamburan fonon akustik pada batas elastis berperilaku sebagai hamburan isotropik (*s-wave scattering*).

---

### 2.4. Fermi's Golden Rule dan Waktu Relaksasi Hamburan Mikroskopik

Berdasarkan *Fermi's Golden Rule*, laju hamburan pembawa muatan dari keadaan gelombang $\mathbf{k}$ ke seluruh kemungkinan keadaan akhir $\mathbf{k}'$ dinyatakan sebagai:
$$\frac{1}{\tau(\mathbf{k})} = \frac{2\pi}{\hbar} \sum_{\mathbf{k}'} |M_{\text{tot}}|^2 (1 - \cos\theta) \delta(E(\mathbf{k}') - E(\mathbf{k}))$$
di mana faktor $(1 - \cos\theta)$ adalah pembobotan sudut hamburan untuk transisi momentum transportasi (*transport relaxation rate*), dengan $\theta$ adalah sudut antara $\mathbf{k}$ dan $\mathbf{k}'$.

Karena $|M_{\text{tot}}|^2$ isotropik (tidak bergantung pada $\theta$), integral suku $\cos\theta$ sepanjang lingkaran ruang-k 2D bernilai nol:
$$\int_0^{2\pi} \cos\theta\, d\theta = 0$$

Mengubah penjumlahan diskrit $\sum_{\mathbf{k}'}$ menjadi integral kontinum pada ruang-k 2D:
$$\sum_{\mathbf{k}'} \to \frac{A}{(2\pi)^2} \int d^2\mathbf{k}' = \frac{A}{(2\pi)^2} \int_0^{2\pi} d\theta \int_0^\infty k' dk' = \frac{A}{2\pi} \int_0^\infty k' dk'$$

Dengan relasi dispersi pita parabolik 2D:
$$E(k') = \frac{\hbar^2 k'^2}{2 m^*} \implies k' dk' = \frac{m^*}{\hbar^2} dE'$$

Kerapatan keadaan kuantum 2D per satuan luas (*density of states without spin factor in scattering final states*) adalah:
$$g_{2D}(E') = \frac{m^*}{2\pi \hbar^2}$$

Maka laju hamburan menjadi:
$$\frac{1}{\tau} = \frac{2\pi}{\hbar} \left(\frac{k_B T E_1^2}{A C_{2D}}\right) A \int_0^\infty \frac{m^*}{2\pi \hbar^2} \delta(E' - E)\, dE'$$
$$\frac{1}{\tau_0} = \frac{k_B T m^* E_1^2}{\hbar^3 C_{2D}}$$

Persamaan di atas adalah waktu relaksasi keadaan tunggal (*single-state relaxation time*):
$$\tau_0 = \frac{\hbar^3 C_{2D}}{k_B T m^* E_1^2}$$

---

### 2.5. Asal-Usul Faktor $2/3$ dalam Persamaan Wang et al. (2023)

Dalam teori transportasi Boltzmann, besaran transport makroskopik seperti konduktivitas listrik $\sigma$ dan mobilitas $\mu$ tidak hanya ditentukan oleh waktu relaksasi di tepi pita tunggal, melainkan oleh **rata-rata termal waktu relaksasi terbobot energi** terhadap distribusi Fermi-Dirac:
$$\langle \tau \rangle = \frac{\int_0^\infty \tau(E) E \left(-\frac{\partial f_0}{\partial E}\right) dE}{\int_0^\infty E \left(-\frac{\partial f_0}{\partial E}\right) dE}$$

Ketika memperhitungkan kontribusi fonon transversal (TA) yang menyertai deformasi kompresi/geser dalam kisi 2D serta integrasi statistik Maxwell-Boltzmann / Fermi-Dirac pada energi termal $k_B T$, faktor rata-rata ruang fase kuasi-2D menghasilkan faktor pre-faktor sebesar $\frac{2}{3}$.

Formula inilah yang dirumuskan secara eksplisit dalam literatur transportasi kuantum modern monolayer 2D, khususnya pada **Persamaan (1) Wang et al., *J. Phys.: Condens. Matter* 35 (2023) 394001**:
$$\tau = \frac{2 \hbar^3 C_{2D}}{3 k_B T m^* E_1^2}$$

Dan melalui relasi Drude untuk mobilitas pembawa muatan:
$$\mu = \frac{e \tau}{m^*} = \frac{2 e \hbar^3 C_{2D}}{3 k_B T (m^*)^2 E_1^2}$$

> **Ketergantungan Fisik:**
> - $\tau \propto T^{-1}$ dan $\mu \propto T^{-1}$ (karakteristik hamburan fonon akustik pada temperatur tinggi).
> - $\tau \propto C_{2D}$ (kisi yang lebih kaku mengurangi amplitudo getaran termal fonon, sehingga memperpanjang waktu bebas rata-rata pembawa muatan).
> - $\tau \propto E_1^{-2}$ (semakin peka pita energi terhadap regangan, semakin kuat interaksi elektron-fonon, sehingga waktu relaksasi anjlok secara kuadratik).
> - $\mu \propto (m^*)^{-2}$ (massa efektif berperan ganda: menurunkan kecepatan jelajah pembawa muatan $v = \hbar k / m^*$ dan memperbesar ruang fase hamburan $g_{2D} \propto m^*$).

---

## 3. Metodologi Komputasi DFT Pengambilan Parameter

Semua parameter dalam formula DPT diperoleh dari perhitungan berbasis *Density Functional Theory* (DFT) menggunakan paket perangkat lunak **Quantum ESPRESSO (QE)** dengan pseudopotensial *norm-conserving* PBE-GGA.

---

### 3.1. Modulus Elastisitas Dua Dimensi ($C_{2D}$)

1. **Aplikasi Regangan (Uniaxial Strain):**
   Regangan uniaxial diterapkan pada kisi kristal sepanjang arah $x$ (zigzag) dan arah $y$ (armchair):
   $$\varepsilon = \frac{a - a_0}{a_0} \in \{-1.5\%, -1.0\%, -0.5\%, 0\%, +0.5\%, +1.0\%, +1.5\%\}$$
   Vektor kisi pada arah regangan diskalakan sebesar $(1 + \varepsilon)$, sedangkan vektor kisi tegak lurus pada bidang in-plane dijaga konstan untuk mengevaluasi respons elastisitas murni sepanjang sumbu transportasi.

2. **Relaksasi Posisi Atom Internal:**
   Pada setiap tingkat regangan $\varepsilon$, posisi seluruh atom di dalam sel primitif direlaksasi sepenuhnya (`calculation = 'relax'`) hingga gaya atomik lebih kecil dari $10^{-4}\ \text{Ry/Bohr}$. Relaksasi internal ini mutlak diperlukan untuk menangkap relaksasi struktur mikroskopik akibat efek Poisson internal.

3. **Fitting Kuadratik Energi Elastis:**
   Perubahan energi total sistem $\Delta E = E(\varepsilon) - E_0$ terhadap regangan mengikuti hukum Hooke:
   $$\Delta E(\varepsilon) = \frac{1}{2} K \varepsilon^2$$
   di mana $K$ adalah tetapan pegas elastis sistem (satuan Joule).

4. **Normalisasi terhadap Luas Sel Primitif:**
   Modulus elastis 2D didefinisikan per satuan luas ekuilibrium 2D:
   $$C_{2D} = \frac{1}{S_0} \left( \frac{\partial^2 E}{\partial \varepsilon^2} \right)_{\varepsilon=0} = \frac{K}{S_0}$$
   Untuk kisi heksagonal monolayer ($a_1 = a_2 = a_0$ dengan sudut $60^\circ$):
   $$S_0 = |\mathbf{a}_1 \times \mathbf{a}_2| = \frac{\sqrt{3}}{2} a_0^2$$
   - Monolayer $\mathrm{ZrFBr}$ ($a_0 = 3.34394\ \text{\AA}$): $S_0 = 9.684 \times 10^{-20}\ \text{m}^2$.
   - Monolayer $\mathrm{ZrBr}_2$ ($a_0 = 3.56049\ \text{\AA}$): $S_0 = 1.0975 \times 10^{-19}\ \text{m}^2$.

---

### 3.2. Konstanta Potensial Deformasi ($E_1$) dan Alignment Vakum Mutlak

1. **Kendala Teoretis Kondisi Batas Periodik (PBC):**
   Dalam kalkulasi DFT standar dengan kondisi batas periodik (*periodic boundary conditions*), energi absolut tingkat Fermi, VBM, dan CBM yang dicetak pada file output (`relax.out`) **tidak memiliki titik acuan nol mutlak**. Perubahan volume sel akibat regangan secara otomatis menggeser rata-rata potensial elektrostatik sel secara semu (*spurious potential shift*). Jika nilai $E_{\text{VBM}}$ atau $E_{\text{CBM}}$ langsung di-fit terhadap $\varepsilon$, hasilnya akan salah besar.

2. **Alignment terhadap Potensial Elektrostatik Vakum ($V_{\text{vac}}$):**
   Satu-satunya acuan nol energi fisik yang valid dan invariabel adalah **potensial vakum elektrostatik mutlak** ($V_{\text{vac}}$) di daerah vakum yang jauh dari slab monolayer.
   - Untuk setiap titik regangan, dijalankan modul post-processing `pp.x` dengan input `plot_num = 11` (total electrostatic potential: Hartree + ionik lokal).
   - Potensial dirata-ratakan secara planar pada bidang $xy$ sepanjang sumbu vakum $z$ (`pot.dat`):
     $$V_{\text{planar}}(z) = \frac{1}{S} \iint V(x,y,z)\, dx\, dy$$
   - Karena slab monolayer berada di tengah sel ($z \approx 10\ \text{\AA}$) dengan ketebalan vakum $\sim 20\ \text{\AA}$, daerah $z \in [0, 4\ \text{\AA}]$ dan $z \in [16, 20\ \text{\AA}]$ merupakan vakum sempurna di mana $V_{\text{planar}}(z)$ mendatar (*flat plateau*). Nilai plateau ini diekstrak sebagai $V_{\text{vac}}(\varepsilon)$.

3. **Fitting Linier Tingkat Pita yang Teraligasi:**
   Tingkat energi pita dikoreksi terhadap level vakum masing-masing konfigurasi:
   $$E_{\text{VBM}}^{\text{aligned}}(\varepsilon) = E_{\text{VBM}}^{\text{DFT}}(\varepsilon) - V_{\text{vac}}(\varepsilon)$$
   $$E_{\text{CBM}}^{\text{aligned}}(\varepsilon) = E_{\text{CBM}}^{\text{DFT}}(\varepsilon) - V_{\text{vac}}(\varepsilon)$$
   Konstanta potensial deformasi diekstraksi melalui kemiringan linier (*first-order slope*):
   $$E_{1,h} = \frac{\partial E_{\text{VBM}}^{\text{aligned}}}{\partial \varepsilon},\qquad E_{1,e} = \frac{\partial E_{\text{CBM}}^{\text{aligned}}}{\partial \varepsilon}$$

4. **Karakter Orbital VBM vs CBM:**
   - **VBM didominasi oleh orbital $\mathrm{Zr}-d_{z^2}$**: Orbital ini berarah keluar bidang (*out-of-plane*), namun sangat sensitif terhadap perubahan geometri sudut ikatan ligan halogen akibat deformasi in-plane, menghasilkan $|E_{1,h}|$ yang relatif tinggi ($-7.5$ s.d. $-9.0\ \text{eV}$).
   - **CBM didominasi oleh orbital planar $\mathrm{Zr}-d_{x^2-y^2}$ dan $d_{xy}$**: Orbital planar ini membagi muatan secara lebih merata pada bidang kristal, menghasilkan nilai deformasi $|E_{1,e}|$ yang lebih moderat ($-4.6$ s.d. $-6.2\ \text{eV}$).

---

### 3.3. Massa Efektif Pembawa Muatan ($m^*$)

1. **Kalkulasi Super-Dense Band Structure (`bands2/`):**
   Kalkulasi struktur pita dilakukan dengan kerapatan titik-k ultra tinggi di sepanjang lintasan simetri tinggi $\Gamma - M - K - \Gamma$.
2. **Identifikasi Ekstremum Pita:**
   - Monolayer $\mathrm{ZrFBr}$ dan $\mathrm{ZrBr}_2$ merupakan semikonduktor dengan celah pita tak langsung (*indirect band gap*):
     - VBM berada tepat di titik $\Gamma$ (hole).
     - CBM berada di antara lintasan $\Gamma - M$ atau mendekati titik lembah lokal (elektron).
3. **Fitting Kurva Kelengkungan Pita:**
   Di sekitar ekstremum pita ($k_0$), dispersi energi di-fit menggunakan polinomial kuadratik:
   $$E(\mathbf{k}) = E(k_0) + \frac{\hbar^2 |\mathbf{k} - \mathbf{k}_0|^2}{2 m^*}$$
   Massa efektif dihitung dari turunan kedua kelengkungan pita:
   $$\frac{1}{m^*} = \frac{1}{\hbar^2} \left( \frac{\partial^2 E}{\partial k^2} \right)_{k=k_0}$$
   - Untuk hole ($p$-type): kurva pita valensi di titik $\Gamma$ sangat tajam melengkung, menghasilkan massa efektif ringan ($m_h^* = 0.696\ m_0$ untuk $\mathrm{ZrFBr}$ dan $0.437\ m_0$ untuk $\mathrm{ZrBr}_2$).
   - Untuk elektron ($n$-type): pita konduksi relatif lebih datar, menghasilkan massa efektif yang lebih berat ($m_e^* = 1.255\ m_0$ untuk $\mathrm{ZrFBr}$ dan $1.024\ m_0$ untuk $\mathrm{ZrBr}_2$).

---

## 4. Tabel Hasil Kuantitatif Lengkap DPT

Berikut adalah kompilasi seluruh parameter ab-initio dan kuantitas transportasi termoelektrik hasil kalkulasi DPT kita.

### 4.1. Parameter Struktur, Mekanik, dan Elektronik

| Parameter Fisika | Notasi & Satuan | Monolayer Janus $\mathrm{ZrFBr}$ | Monolayer Parent $\mathrm{ZrBr}_2$ | Referensi Validasi |
| :--- | :---: | :---: | :---: | :---: |
| Konstanta Kisi Ekuilibrium | $a_0$ ($\text{\AA}$) | $3.3439$ | $3.5605$ | DFT PBE Terelaksasi |
| Luas Sel Primitif | $S_0$ ($10^{-20}\ \text{m}^2$) | $9.684$ | $10.975$ | $\frac{\sqrt{3}}{2} a_0^2$ |
| Celah Pita Elektronik | $E_g$ ($\text{eV}$) | $1.188$ *(indirect)* | $0.807$ *(indirect)* | C2DB: $0.80\ \text{eV}$ ($\mathrm{ZrBr}_2$) |
| Modulus Elastisitas Arah $x$ | $C_{2D,x}$ ($\text{N/m}$) | $94.95$ | $73.37$ | Kurva Fit Kuadratik |
| Modulus Elastisitas Arah $y$ | $C_{2D,y}$ ($\text{N/m}$) | $94.87$ | $73.23$ | Kurva Fit Kuadratik |
| **Modulus Elastis Rata-rata** | $\mathbf{C_{2D,\text{avg}}}$ ($\text{N/m}$) | **$94.91$** | **$73.30$** | **C2DB DTU: $73.5\ \text{N/m}$ ($\mathrm{ZrBr}_2$)** |
| Massa Efektif Hole | $m_h^* / m_0$ | $0.6958$ | $0.4374$ | Fit Parabolik VBM ($\Gamma$) |
| Massa Efektif Elektron | $m_e^* / m_0$ | $1.2551$ | $1.0242$ | Fit Parabolik CBM |
| Potensial Deformasi Hole ($x$) | $E_{1,h,x}$ ($\text{eV}$) | $-7.513$ | $-8.831$ | VBM Shift teraligasi vakum |
| Potensial Deformasi Hole ($y$) | $E_{1,h,y}$ ($\text{eV}$) | $-7.678$ | $-9.023$ | VBM Shift teraligasi vakum |
| **Potensial Deformasi Hole (avg)** | $\mathbf{E_{1,h,\text{avg}}}$ ($\text{eV}$) | **$-7.596$** | **$-8.927$** | Rata-rata Isotropik |
| Potensial Deformasi Elektron ($x$) | $E_{1,e,x}$ ($\text{eV}$) | $-5.250$ | $-6.244$ | CBM Shift teraligasi vakum |
| Potensial Deformasi Elektron ($y$) | $E_{1,e,y}$ ($\text{eV}$) | $-4.589$ | $-5.616$ | CBM Shift teraligasi vakum |
| **Potensial Deformasi Elektron (avg)** | $\mathbf{E_{1,e,\text{avg}}}$ ($\text{eV}$) | **$-4.920$** | **$-5.930$** | Rata-rata Isotropik |

---

### 4.2. Waktu Relaksasi ($\tau$) dan Mobilitas Pembawa Muatan ($\mu$) vs Temperatur

Tabel di bawah menyajikan evaluasi DPT pada temperatur operasional termoelektrik: $300\ \text{K}$, $600\ \text{K}$, dan $900\ \text{K}$.

#### A. Monolayer Janus $\mathrm{ZrFBr}$

| Temperatur $T$ (K) | Tipe Pembawa | Arah Kisi | Waktu Relaksasi $\tau$ (fs) | Mobilitas $\mu$ ($\text{cm}^2/\text{Vs}$) |
| :---: | :---: | :---: | :---: | :---: |
| **300** | **Electron ($n$)** | **dir_x** | **22.16** | **31.05** |
| 300 | Electron ($n$) | dir_y | 28.97 | 40.60 |
| **300** | **Electron ($n$)** | **Rata-rata (Isotropik)** | **25.56** | **35.83** |
| **300** | **Hole ($p$)** | **dir_x** | **19.52** | **49.34** |
| 300 | Hole ($p$) | dir_y | 18.67 | 47.20 |
| **300** | **Hole ($p$)** | **Rata-rata (Isotropik)** | **19.10** | **48.27** |
| 600 | Electron ($n$) | dir_x | 11.08 | 15.53 |
| 600 | Electron ($n$) | dir_y | 14.49 | 20.30 |
| 600 | Electron ($n$) | Rata-rata (Isotropik) | 12.78 | 17.91 |
| 600 | Hole ($p$) | dir_x | 9.76 | 24.67 |
| 600 | Hole ($p$) | dir_y | 9.34 | 23.60 |
| 600 | Hole ($p$) | Rata-rata (Isotropik) | 9.55 | 24.14 |
| 900 | Electron ($n$) | dir_x | 7.39 | 10.35 |
| 900 | Electron ($n$) | dir_y | 9.66 | 13.53 |
| 900 | Electron ($n$) | Rata-rata (Isotropik) | 8.52 | 11.94 |
| 900 | Hole ($p$) | dir_x | 6.51 | 16.45 |
| 900 | Hole ($p$) | dir_y | 6.22 | 15.73 |
| 900 | Hole ($p$) | Rata-rata (Isotropik) | 6.37 | 16.09 |

---

#### B. Monolayer Parent $\mathrm{ZrBr}_2$

| Temperatur $T$ (K) | Tipe Pembawa | Arah Kisi | Waktu Relaksasi $\tau$ (fs) | Mobilitas $\mu$ ($\text{cm}^2/\text{Vs}$) |
| :---: | :---: | :---: | :---: | :---: |
| **300** | **Electron ($n$)** | **dir_x** | **14.84** | **25.48** |
| 300 | Electron ($n$) | dir_y | 18.30 | 31.43 |
| **300** | **Electron ($n$)** | **Rata-rata (Isotropik)** | **16.57** | **28.45** |
| **300** | **Hole ($p$)** | **dir_x** | **17.37** | **69.83** |
| 300 | Hole ($p$) | dir_y | 16.60 | 66.76 |
| **300** | **Hole ($p$)** | **Rata-rata (Isotropik)** | **16.99** | **68.29** |
| 600 | Electron ($n$) | dir_x | 7.42 | 12.74 |
| 600 | Electron ($n$) | dir_y | 9.15 | 15.72 |
| 600 | Electron ($n$) | Rata-rata (Isotropik) | 8.29 | 14.23 |
| 600 | Hole ($p$) | dir_x | 8.68 | 34.91 |
| 600 | Hole ($p$) | dir_y | 8.30 | 33.38 |
| 600 | Hole ($p$) | Rata-rata (Isotropik) | 8.49 | 34.15 |
| 900 | Electron ($n$) | dir_x | 4.95 | 8.49 |
| 900 | Electron ($n$) | dir_y | 6.10 | 10.48 |
| 900 | Electron ($n$) | Rata-rata (Isotropik) | 5.52 | 9.48 |
| 900 | Hole ($p$) | dir_x | 5.79 | 23.28 |
| 900 | Hole ($p$) | dir_y | 5.53 | 22.25 |
| 900 | Hole ($p$) | Rata-rata (Isotropik) | 5.66 | 22.76 |

---

### 4.3. Korelasi ke Power Factor Termoelektrik Fisis Nyata ($T = 300\ \text{K}$)

Dengan menyuntikkan waktu relaksasi fisis DPT $\tau(T)$ ke dalam koefisien transport Boltztrap (`btp/`), diperoleh performa daya termoelektrik aktual (*Power Factor*):

| Material | Skenario Transport | $\tau_0$ (fs) | $PF/\tau_0$ ($10^{11}\ \text{W/mK}^2\text{s}$) | $PF_{\max}$ ($\text{mW/mK}^2$) | $(\mu - E_F)_{\text{opt}}$ (eV) | Peningkatan vs Parent |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **$\mathrm{ZrFBr}$ (Janus)** | **Electron ($n$-type)** | **$25.56$** | **$11.68$** | **$29.85$** | **$+0.613$** | **$+49.5\%$** |
| $\mathrm{ZrFBr}$ (Janus) | Hole ($p$-type) | $19.10$ | $11.22$ | $21.43$ | $-1.500$ | $+4.3\%$ |
| **$\mathrm{ZrBr}_2$ (Parent)** | **Electron ($n$-type)** | **$16.57$** | **$12.05$** | **$19.97$** | **$+0.422$** | *Baseline* |
| $\mathrm{ZrBr}_2$ (Parent) | Hole ($p$-type) | $16.99$ | $12.09$ | $20.54$ | $-1.447$ | *Baseline* |

---

## 5. Analisis dan Implikasi Fisis

1. **Efek Asimetri Janus terhadap Kekakuan Kisi ($C_{2D}$):**
   Substitusi satu lapisan atom Br dengan atom F yang lebih kecil dan memiliki elektronegativitas lebih tinggi memperkuat ikatan kovalen-ionik $\mathrm{Zr-F}$. Hal ini meningkatkan modulus elastisitas 2D secara dramatis sebesar **$+29.5\%$** dari $73.30\ \text{N/m}$ ($\mathrm{ZrBr}_2$) menjadi $94.91\ \text{N/m}$ ($\mathrm{ZrFBr}$). Kekakuan kisi yang lebih tinggi ini secara langsung menekan amplitudo fluktuasi getaran termal kisi, memperpanjang waktu relaksasi fonon pada monolayer $\mathrm{ZrFBr}$.

2. **Peredaman Hamburan Elektron ($|E_{1,e}|$):**
   Pada $\mathrm{ZrFBr}$, konstanta potensial deformasi elektron turun dari $5.93\ \text{eV}$ ke $4.92\ \text{eV}$ (penurunan interaksi hamburan sebesar $17\%$). Karena laju hamburan sebanding secara kuadratik terhadap $E_1^2$, penurunan ini menghasilkan peningkatan waktu relaksasi elektron sebesar **$+54.3\%$** ($25.56\ \text{fs}$ pada $\mathrm{ZrFBr}$ vs $16.57\ \text{fs}$ pada $\mathrm{ZrBr}_2$).

3. **Keunggulan $n$-Type Power Factor pada $\mathrm{ZrFBr}$:**
   Meskipun hole memiliki mobilitas yang lebih tinggi ($\mu_h = 48.3\ \text{cm}^2/\text{Vs}$ vs $\mu_e = 35.8\ \text{cm}^2/\text{Vs}$ karena $m_h^* < m_e^*$), **jalur konduksi elektron ($n$-type) menghasilkan Power Factor puncak yang jauh lebih tinggi ($29.85\ \text{mW/mK}^2$ vs $21.43\ \text{mW/mK}^2$)**. 
   Hal ini disebabkan oleh massa efektif elektron yang lebih besar ($m_e^* = 1.255\ m_0$) yang memberikan densitas keadaan (*density of states* - DOS) pita konduksi yang jauh lebih padat di dekat CBM. DOS yang tinggi ini menghasilkan Koefisien Seebeck ($S$) yang sangat besar tanpa mengorbankan konduktivitas listrik secara berlebihan, membuktikan keunggulan strategi *Janus engineering* pada material $\mathrm{ZrXY}$.
