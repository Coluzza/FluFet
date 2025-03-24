#!/bin/bash
#set -x
PROT1="$1"

#Parse protein name and replace - with _
GO_PROT1=`basename ${PROT1%.*} -out| awk '{gsub("-","_",$0); print $0}'`

#define variable to read data from
READDATA=`basename ${PROT1%.*}`

#insert head of the data file
cat > in.${GO_PROT1}_graft  <<\EOF

units		lj
atom_style	full
atom_modify	first forcegroup #needed for SRDmod
dimension	3
boundary	p p p
EOF


echo "read_data ${GO_PROT1}-out.atomparam" >> in.${GO_PROT1}_graft

cat >> in.${GO_PROT1}_graft  <<\EOF
group dummy type 1 2
group protein subtract all dummy

#Shift polymer to the middle of the box, since box size change will rip bonds apart otherwise
#Always choose the box boundaries of the input xlo 0 ylo 0 zlo 0 otherwise you might get problems when the old box is not included in the new one
variable xshift1 equal -xcm(protein,x)+(xhi-xlo)/2
variable yshift1 equal -xcm(protein,y)+(yhi-ylo)/2
variable zshift1 equal -xcm(protein,z)+(zhi-zlo)/2
displace_atoms protein move ${xshift1} ${yshift1} ${zshift1}

#Define the geometry you want
variable        xlo equal 0.0
variable        xhi equal 20
variable        ylo equal 0.0
variable        yhi equal 20
variable        zlo equal 0.0
variable        zhi equal 20

#Change box to that geometry
change_box all x final ${xlo} ${xhi} y final ${ylo} ${yhi} z final ${zlo} ${zhi} boundary p p p 

#Shift protein to the center of the cylindrical region again
variable xshift2 equal -xcm(protein,x)+(xhi-xlo)/2
variable yshift2 equal -xcm(protein,y)+(yhi-ylo)/2
variable zshift2 equal -xcm(protein,z)+(zhi-zlo)/2
displace_atoms protein move ${xshift2} ${yshift2} ${zshift2}
set group protein image 0 0 0

#Spacing for flow profile if needed
variable 		xspacing equal 10.0

#Geometry of the cylinder
variable        middleyz equal ${yhi}/2.0
variable        radius equal 7
variable        inradius equal  ${radius}-0.5
#dummy variables to have no wall in x-direction
variable 	xlo2 equal ${xlo}-4.0
variable 	xhi2 equal ${xhi}+4.0

#Create regions of the cylinder
region          poreinner cylinder x ${middleyz} ${middleyz} ${inradius} ${xlo} ${xhi}
region 		porewall cylinder x ${middleyz} ${middleyz} ${radius} ${xlo2} ${xhi2}


#Define variables of the go model
variable alpha equal 2.0*3.8
variable sigma equal 1.0

variable eps_h equal 4.0
variable bonds_NN equal ${sigma}*1.12246204830938
variable cutoff equal ${bonds_NN}*4.0

variable GO equal -4
variable GO2 equal 4

variable prot_k equal ${eps_h}*50 #changed this to 200 from 20
variable angle_k equal ${eps_h}*6
variable angle_k2 equal ${eps_h}*6
variable torsional_1 equal ${eps_h}*6
variable torsional_2 equal ${eps_h}*5.4

variable cutlj equal 1.12

neighbor	    0.3 bin
neigh_modify	delay 0 every 1 check yes

EOF

#Insert bonds and angles and dihedrals
cat ${GO_PROT1}-out.int_rescale >> in.${GO_PROT1}_graft


cat >> in.${GO_PROT1}_graft  <<\EOF

#create atoms at the wall
create_atoms 2 single 10.0 10.0 3.0

group wallpoint type 2
group forcegroup union protein wallpoint

#use a soft potential to pull residue of choice to the wall point
pair_style       hybrid/overlay lj/cut ${bonds_NN} sigmoidal ${cutoff} soft 10.0
pair_coeff * * lj/cut 1 ${sigma}
pair_coeff 1 * lj/cut 0 0 0.0
pair_coeff * * soft 0 0.0
pair_coeff 2 3 soft -100 10.0
pair_coeff * * sigmoidal 0 0 0 0
pair_modify shift yes
special_bonds	lj 0.0 1.0 1.0


EOF



cat ${GO_PROT1}-out.paircoeff >> in.${GO_PROT1}_graft


cat >> in.${GO_PROT1}_graft  <<\EOF

fix e1 forcegroup nve
fix e2 forcegroup langevin 1.0 1.0 1 12231 zero yes
fix f1 wallpoint setforce 0.0 0.0 0.0
fix		w1 protein wall/region porewall lj126 1.0 1.0 ${cutlj} #pushes protein away from the wall


dump d1 protein dcd 100 trajprotein.dcd
dump d2 wallpoint atom 100 wallpoint.lammpstrj

timestep 0.002
thermo 1000
run 10000

#create fixed bond
group ankerres type 3
create_bonds many ankerres wallpoint 1 0.5 2.0

run 10000

unfix e1
unfix e2
reset_timestep 0

EOF



#cat ${GO_PROT1}-out.int_rescale >> in.${GO_PROT1}_graft


cat >> in.${GO_PROT1}_graft  <<\EOF

pair_style       hybrid/overlay lj/cut ${bonds_NN} sigmoidal ${cutoff}
pair_coeff * * lj/cut 1 ${sigma}
pair_coeff 1 * lj/cut 0 0 0.0
pair_coeff * * sigmoidal 0 0 0 0
pair_modify shift yes
special_bonds	lj 0.0 1.0 1.0


EOF



cat ${GO_PROT1}-out.paircoeff >> in.${GO_PROT1}_graft



#insert end commands
cat >> in.${GO_PROT1}_graft  <<\EOF

#create srd solvent
lattice        sc 1.0
lattice        sc 10.0
create_atoms   1 region poreinner

group solvent type 1
group moving union solvent protein

#initial velocity of protein can be changed
velocity	solvent create 1.0 44754895 
velocity        protein set NULL NULL NULL units box

neigh_modify	include forcegroup
comm_modify	mode multi group forcegroup vel yes


timestep	0.002

fix             s1 moving srd 5 NULL 1.0 1.0 49894 rescale no tstat yes shift yes 1234 beadgroup protein  fflow 0.1 0.1 rcyl ${inradius}  virtual 10 alpha 120


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

run		1000000

EOF
