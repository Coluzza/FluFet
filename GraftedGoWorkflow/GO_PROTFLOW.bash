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

neighbor	    0.3 multi
neigh_modify	delay 0 every 10 check yes
neigh_modify	include forcegroup
comm_modify	mode multi cutoff/multi * 80 group forcegroup vel yes


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
cat ${GO_PROT1}-out.pullforce >> in.${GO_PROT1}_lmp

#fix dock anchor addforce 0.0 0.0 -60.0
#fix dock2 min addforce 0.0 0.0 -50.0
cat >> in.${GO_PROT1}_lmp  <<\EOF
dump 1 all custom  1000  4_docking.lammpstrj  id  type  x y z  ix iy iz
fix 2 all nve
fix 3 all langevin 1 1 1 ${seed} zero yes
#fix 5 all momentum 1 linear 1 1 1 angular rescale
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

cat >> in.${GO_PROT1}_lmp  <<\EOF





EOF

#exit


echo 'variable prot_k equal ${eps_h}*20'  >> in.${GO_PROT1}_lmp
echo 'variable linker_k equal ${eps_h}*20'  >> in.${GO_PROT1}_lmp
cat ${GO_PROT1}-out.harmonicbonds >> in.${GO_PROT1}_lmp

cat ${GO_PROT1}-out.srd >> in.${GO_PROT1}_lmp


#insert end commands
cat >> in.${GO_PROT1}_lmp  <<\EOF
group moving region poreinner 
#initial velocity of protein can be changed
velocity	solvent create ${ref_temp} ${seed} 
velocity  protein set NULL NULL NULL units box




unfix 2 
unfix 3 
timestep	0.002

fix             s1 moving srd 5 NULL 1.0 1.0 49894 rescale no tstat yes shift yes 1234 beadgroup free  fflow 0.0001 0.0001 rcyl ${inradius}  virtual 10 alpha 120


#just for orientation
compute cm1 protein com
compute rg protein gyration

#control variables
thermo_style custom step temp ebond eangle edihed epair pe ke etotal c_cm1[1] c_rg
thermo		1000

#write out gyration tensor
fix             r1 protein ave/time 1000 1 1000 c_rg[1] c_rg[2] c_rg[3] c_rg[4] c_rg[5] c_rg[6] file tensor.data

#write out flow profile if needed
#reset_timestep 0
#compute 		cc1 solvent chunk/atom bin/3d x 0.0 ${xspacing} y 0.0 0.2 z 0.0 0.2 units box
#fix             	v1 solvent ave/chunk 100 10000 1000000 cc1 vx vy vz temp norm all file vel_srd_f1t

#dump d2 solvent dcd 1000 solv.dcd
#restart 50000000 re1.data re2.data 

dump 1 forcegroup custom  1000  7_flow.lammpstrj  id  type  x y z  ix iy iz
#run		1000000

EOF
