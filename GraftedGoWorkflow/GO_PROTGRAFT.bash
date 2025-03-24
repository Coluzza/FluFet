#!/bin/bash
#set -x
PROT1="$1"

#Parse protein name and replace - with _
GO_PROT1=`basename ${PROT1%.*} -out| awk '{gsub("-","_",$0); print $0}'`

#define variable to read data from
READDATA=`basename ${PROT1%.*}`

#insert head of the data file
cat > in.${GO_PROT1}_lmp  <<\EOF

units		lj
atom_style	full
atom_modify	first forcegroup #needed for SRDmod
dimension	3
boundary	p p p
EOF


echo "read_data ${GO_PROT1}-out.atomparam" >> in.${GO_PROT1}_lmp

cat ${GO_PROT1}-out.groups >> in.${GO_PROT1}_lmp

cat ${GO_PROT1}-out.variables >> in.${GO_PROT1}_lmp

#Insert bonds and angles and dihedrals
cat ${GO_PROT1}-out.harmonicbonds >> in.${GO_PROT1}_lmp
cat ${GO_PROT1}-out.GOdihedral >> in.${GO_PROT1}_lmp



cat >> in.${GO_PROT1}_lmp  <<\EOF

neighbor	    0.3 bin
neigh_modify	delay 0 every 10 check yes
neigh_modify	include forcegroup
comm_modify	mode single cutoff 80 group forcegroup vel yes



EOF

cat ${GO_PROT1}-out.fix >> in.${GO_PROT1}_lmp

cat >> in.${GO_PROT1}_lmp  <<\EOF



velocity      all create 1  ${seed} 
pair_style soft ${cutoff_C}
pair_coeff * * 10.0
variable prefactor equal ramp(10,100)
pair_modify shift yes
fix 4 all adapt 1 pair soft a * * v_prefactor

fix 2 all nve
fix 3 all langevin 1 1 1 ${seed} zero yes
fix 5 all momentum 1 linear 1 1 1 angular rescale
dump 1 all custom 10 1_minimize.lammpstrj id  type  x y z  ix iy iz
thermo          100
timestep        0.002
run 2000
undump 1
unfix 4
dump 1 all custom 1000 2_unfolding.lammpstrj id  type  x y z  ix iy iz
velocity      all create 10.0  ${seed} 

pair_style lj/sf ${cutoff_C}
pair_coeff * * 1 ${sigma_C}
pair_modify shift yes

fix 2 all nve
fix 3 all langevin 10.0 10.0 1 ${seed} zero yes
fix 5 all momentum 10 linear 1 1 1 angular rescale
thermo          100
timestep        0.002
run 2000
undump 1

velocity      all create 1  ${seed} 
reset_timestep 0
EOF


echo 'variable linker_k equal ${eps_h}*5'  >> in.${GO_PROT1}_lmp
cat ${GO_PROT1}-out.harmonicbonds >> in.${GO_PROT1}_lmp
cat ${GO_PROT1}-out.intparam >> in.${GO_PROT1}_lmp


cat >> in.${GO_PROT1}_lmp  <<\EOF
compute cc1 all pair sigmoidal evdwl

thermo_style custom step temp ebond eangle edihed epair pe ke etotal c_cc1
thermo 1000

dump 1 all custom  100  3_folding.lammpstrj  id  type  x y z  ix iy iz
fix 2 all nve
fix 3 all langevin 1 1 1 ${seed} zero yes
fix 5 all momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 10000
undump 1
EOF
#cat ${GO_PROT1}-out.pullforcea >> in.${GO_PROT1}_lmp


cat >> in.${GO_PROT1}_lmp  <<\EOF
fix dock anchor addforce 0.0 0.0 -50.0
fix dock2 min addforce 0.0 0.0 -50.0
dump 1 all custom  1000  4a_docking.lammpstrj  id  type  x y z  ix iy iz
fix 2 all nve
fix 3 all langevin 1 1 1 ${seed} zero yes
fix 5 all momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 200000
undump 1

fix dock anchor addforce 0.0 0.0 -15.0
fix dock2 min addforce 0.0 0.0 -10.0

dump 1 all custom  1000  4b_docking.lammpstrj  id  type  x y z  ix iy iz
fix 2 all nve
fix 3 all langevin 1 1 1 ${seed} zero yes
fix 5 all momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 200000
undump 1
unfix dock2
EOF


cat ${GO_PROT1}-out.unfix >> in.${GO_PROT1}_lmp
cat >> in.${GO_PROT1}_lmp  <<\EOF
fix 2 free nve
fix 3 free langevin 1 1 1 12231 zero yes
#fix 5 min momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 200000
undump 1
EOF


cat ${GO_PROT1}-out.fixbrooks >> in.${GO_PROT1}_lmp
cat >> in.${GO_PROT1}_lmp  <<\EOF
unfix 5
dump 1 all custom  1000  6_bound_lowT.lammpstrj  id  type  x y z  ix iy iz
fix 2 free nve
fix 3 free langevin 1 1 1 ${seed} zero yes
#fix 5 min momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 1000000
undump 1

dump 1 all custom  1000  6_bound_highT.lammpstrj  id  type  x y z  ix iy iz
fix 2 free nve
fix 3 free langevin ${ref_temp} ${ref_temp} 1 ${seed} zero yes
#fix 5 min momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 1000000
undump 1
EOF
