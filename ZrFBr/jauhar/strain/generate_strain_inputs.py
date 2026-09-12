import os
import numpy as np

# Basis parameter ZrFBr
a0 = 3.343942779 # Angstrom
c0 = 10.000000000 # Angstrom

strains = [-0.015, -0.010, -0.005, 0.000, 0.010, 0.015, 0.020]
# Gunakan 7 titik simetris: -1.5%, -1.0%, -0.5%, 0.0%, +0.5%, +1.0%, +1.5%
strains = [-0.015, -0.010, -0.005, 0.000, 0.005, 0.010, 0.015]

directions = ['x', 'y'] # Armchair (x) dan Zigzag (y)

template_scf = """&CONTROL
  calculation  = 'relax',
  pseudo_dir   = '../../../pseudo/',
  outdir       = './tmp',
  prefix       = 'ZrFBr',
  forc_conv_thr = 1.0d-4,
  etot_conv_thr = 1.0d-5,
/

&SYSTEM
  ibrav       = 0,
  nat         = 3,
  ntyp        = 3,
  ecutwfc     = 65.0,
  nbnd        = 36,
/

&ELECTRONS
  conv_thr    = 1.0d-8,
/

&IONS
/

ATOMIC_SPECIES
  Zr  91.224  Zr.pbe-spn-kjpaw_psl.1.0.0.UPF
  F   18.998  F.pbe-n-kjpaw_psl.1.0.0.UPF
  Br  79.904  Br.pbe-n-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS (angstrom)
{v1x:14.9f} {v1y:14.9f} {v1z:14.9f}
{v2x:14.9f} {v2y:14.9f} {v2z:14.9f}
{v3x:14.9f} {v3y:14.9f} {v3z:14.9f}

ATOMIC_POSITIONS (angstrom)
  Zr   1.671971388   0.965313131  10.197570978
  F    3.343942780   1.930626264  11.502363216
  Br   3.343942780   1.930626264   8.300065808

K_POINTS (automatic)
  16 16 1 0 0 0
"""

template_pp = """&INPUTPP
  prefix   = 'ZrFBr',
  outdir   = './tmp',
  filplot  = 'pot.dat',
  plot_num = 11,
/
&PLOT
  nfile         = 1,
  filepp(1)     = 'pot.dat',
  weight(1)     = 1.0,
  iflag         = 1,
  output_format = 0,
  e1(1) = 0.0, e1(2) = 0.0, e1(3) = 1.0,
  x0(1) = 0.0, x0(2) = 0.0, x0(3) = 0.0,
  nx            = 500,
/
"""

for d in directions:
    dir_path = f"dir_{d}"
    os.makedirs(dir_path, exist_ok=True)
    
    for eps in strains:
        tag = f"eps_{eps:+.3f}"
        sub_dir = os.path.join(dir_path, tag)
        os.makedirs(sub_dir, exist_ok=True)
        
        # Hitung vektor kisi ter-deformasi
        if d == 'x': # Uniaxial strain sepanjang x (armchair)
            v1x = a0 * (1.0 + eps)
            v1y = 0.0
            v1z = 0.0
            
            v2x = (a0 * (1.0 + eps)) / 2.0
            v2y = (a0 * np.sqrt(3.0)) / 2.0
            v2z = 0.0
        elif d == 'y': # Uniaxial strain sepanjang y (zigzag)
            v1x = a0
            v1y = 0.0
            v1z = 0.0
            
            v2x = a0 / 2.0
            v2y = ((a0 * np.sqrt(3.0)) / 2.0) * (1.0 + eps)
            v2z = 0.0
            
        v3x = 0.0
        v3y = 0.0
        v3z = c0
        
        scf_content = template_scf.format(
            v1x=v1x, v1y=v1y, v1z=v1z,
            v2x=v2x, v2y=v2y, v2z=v2z,
            v3x=v3x, v3y=v3y, v3z=v3z
        )
        
        with open(os.path.join(sub_dir, "relax.in"), "w") as f:
            f.write(scf_content)
            
        with open(os.path.join(sub_dir, "pp.in"), "w") as f:
            f.write(template_pp)

print("Selesai men-generate seluruh file input regangan untuk ZrFBr (arah x dan arah y)!")
