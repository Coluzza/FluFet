#!/bin/bash
#set -x
INPUT="$1"

#Parse protein name and replace - with _
PREFIX=`basename ${INPUT%.*} -out| awk '{gsub("-","_",$0); print $0}'`

#define variable to read data from
READDATA=`basename ${INPUT%.*}`

#insert head of the data file
cat > in.${PREFIX}_lmp  <<\EOF

units		lj
atom_style	full
atom_modify	first forcegroup #needed for SRDmod
dimension	3
boundary	p p p
EOF


echo "read_data ${PREFIX}.atomparam" >> in.${PREFIX}_lmp

cat ${PREFIX}.groups >> in.${PREFIX}_lmp

cat ${PREFIX}.variables >> in.${PREFIX}_lmp

#Insert bonds and angles and dihedrals

cat ${PREFIX}.bondcoeff >> in.${PREFIX}_lmp

cat >> in.${PREFIX}_lmp  <<\EOF

neighbor	    0.5 bin
neigh_modify	delay 0 every 10 check yes
neigh_modify	include forcegroup
comm_modify	mode multi cutoff/multi * 50 group forcegroup vel yes



EOF

cat ${PREFIX}.fix >> in.${PREFIX}_lmp

cat >> in.${PREFIX}_lmp  <<\EOF


fix rigid1 all rigid group 1 nanoparticle
neigh_modify exclude group nanoparticle nanoparticle
neigh_modify exclude group anchor_surf anchor_surf

velocity      forcegroup create 1.0 ${seed}
fix rigid2 anchor_surf setforce 0 0 0
velocity  anchor_surf set 0 0 0 units box


################### SOFT REPULSION TO EQUILIBRATE##################
pair_style soft ${cutoff_polymer}
pair_coeff * * 10.0
variable prefactor equal ramp(10,100)
pair_modify shift yes
fix 4 forcegroup adapt 1 pair soft a * * v_prefactor

fix 2 all nve
fix 3 forcegroup langevin 1.0 1.0 1 ${seed} zero yes
fix 5 forcegroup momentum 1 linear 1 1 1 angular rescale
dump 1 output custom 10 1_soft.lammpstrj id  type  x y z  ix iy iz
thermo          1
timestep        0.002
run 2000
undump 1
unfix 4

################### HARD REPULSIONS TO EQUILIBRATE##################

dump 1 output custom 1000 2_hard.lammpstrj id  type  x y z  ix iy iz

pair_style lj/sf ${cutoff_polymer}
pair_coeff * * 1 ${sigma_polymer}
pair_modify shift yes

fix 2 all nve
fix 3 forcegroup langevin 1.0 1.0 1 ${seed} zero yes
fix 5 forcegroup momentum 10 linear 1 1 1 angular rescale
thermo          100
timestep        0.002
run 2000
undump 1

reset_timestep 0
################### REAL INTERACTIONS TO EQUILIBRATE##################
EOF

cat ${PREFIX}.intparam >> in.${PREFIX}_lmp


cat >> in.${PREFIX}_lmp  <<\EOF

thermo_style custom step temp ebond epair pe ke etotal 
thermo 1000

dump 1 output custom  100  3_equil.lammpstrj  id  type  x y z  ix iy iz
fix 2 all nve
fix 3 forcegroup langevin 1.0 1.0 1 ${seed} zero yes
fix 5 forcegroup momentum 1 linear 1 1 1 angular rescale
thermo          1000
timestep        0.002
run 200000
undump 1
unfix 2 
unfix 3 
unfix 5

################### FLOW ##################
EOF


cat ${PREFIX}.srd >> in.${PREFIX}_lmp


#insert end commands
cat >> in.${PREFIX}_lmp  <<\EOF
#initial velocity of protein can be changed
velocity	solvent create 1.0 ${seed} 
#velocity  forcegroup set NULL NULL NULL units box

timestep	0.002

fix             s1 moving srd 5 NULL 1.0 1.0 49894 rescale yes tstat yes shift yes 1234 beadgroup flowgroup  fflow 0.01 0.01 rcyl ${radius} virtual 10 alpha 120
#control variables
#thermo_style custom step temp ebond eangle edihed epair pe ke etotal c_cm1[1] c_rg
thermo		1000


#write out flow profile if needed
reset_timestep 0
compute 		cc1 solvent chunk/atom bin/3d x 0.0 5 y 0.0 0.2 z 0.0 0.2 units box
fix         v1 solvent ave/chunk 10 10000 100000 cc1 vx vy vz temp norm all file vel_srd_f1t

#dump d2 solvent dcd 1000 solv.dcd
#restart 50000000 re1.data re2.data 

dump 1000 output custom  1000  4_flow.lammpstrj  id  type  x y z  ix iy iz
run		1000000

EOF
