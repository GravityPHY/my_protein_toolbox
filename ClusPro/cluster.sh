/usr4/spclpgm/imhaoyu/.local/bin/sblu pdb prep --smod R --out-prefix rec rec.pdb
/usr4/spclpgm/imhaoyu/.local/bin/sblu pdb prep --smod L --out-prefix lig lig.pdb
coeffs=('000' '002' '004' '006')
for coeff in ${coeffs[@]}
do
  /usr4/spclpgm/imhaoyu/.local/bin/sblu measure pwrmsd -n 1500 --only-CA --only-interface \
  --rec rec_nmin.pdb -o clustermat.${coeff}.00 lig_nmin.pdb \
  ft.${coeff}.00 /projectnb2/docking/imhaoyu/ClusPro_TM/rot_mat/rot70k.0.0.6.jm.max_10deg_z-axis_tilt.mol2

  /usr4/spclpgm/imhaoyu/.local/bin/sblu docking cluster -o clustermat.${coeff}.00.clusters \
  clustermat.${coeff}.00

# minimize

  /usr4/spclpgm/imhaoyu/.local/bin/sblu cluspro minimize rec_nmin.pdb rec.psf lig_nmin.pdb lig.psf \
  clustermat.${coeff}.00.clusters ft.${coeff}.00 \
  /projectnb2/docking/imhaoyu/ClusPro_TM/rot_mat/rot70k.0.0.6.jm.max_10deg_z-axis_tilt.mol2 -o model.${coeff}
done