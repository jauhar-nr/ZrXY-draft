#!/bin/bash
# Script untuk menjalankan BoltzTraP2 murni secara lokal
# Menggunakan data dari folder dos/tmp/ yang sudah ada

echo "Memulai interpolasi BoltzTraP2 untuk ZrBr2..."
# -m 8: Multiplier 8
# -e -0.35 -E 0.35: Rentang "Set & Forget" super aman untuk semua jenis material (~10 eV)
btp2 -vv interpolate -m 8 -e -0.35 -E 0.35 ../dos/tmp/

echo "---------------------------------------------------"
echo "Memulai integrasi BoltzTraP2..."
# Integrasi rentang suhu 300K sampai 900K dengan step 300K
btp2 -vv integrate interpolation.bt2 -b 10000 300:1000:300

echo "Perhitungan BoltzTraP2 Selesai! Data disimpan di file interpolation.bt2"
