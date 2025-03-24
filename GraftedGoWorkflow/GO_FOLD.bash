#!/bin/bash
#set -x
PROT1="$1"

#Parse protein name and replace - with _
GO_PROT1=`basename ${PROT1%.*} -out| awk '{gsub("-","_",$0); print $0}'`

#define variable to read data from
READDATA=`basename ${PROT1%.*}`

#insert head of the data file
cat > in.${GO_PROT1}_lmp  <<\EOF

units           lj
atom_style              full
boundary                p p p

variable alpha equal 2.0*3.8
variable sigma equal 1.0
EOF

echo "read_data ${GO_PROT1}-out.atom_rescale" >> in.${GO_PROT1}_lmp 


cat >> in.${GO_PROT1}_lmp  <<\EOF
variable eps_h equal 4.0
variable bonds_NN equal ${sigma}*1.12246204830938
variable cutoff equal ${bonds_NN}*4.0

variable GO equal -4
variable GO2 equal 4

variable prot_k equal ${eps_h}*5
variable angle_k equal ${eps_h}*6
variable angle_k2 equal ${eps_h}*6
variable torsional_1 equal ${eps_h}*6
variable torsional_2 equal ${eps_h}*5.4

variable mytau equal sqrt(1/${eps_h})
variable mydamp equal 1/${mytau}
#variable mydt equal 0.025*
variable mydt equal 0.002*${mytau}
variable cutlj equal 1.12
print "mydamp, mydt ${mydamp} ${mydt}"

#neighbor                 0.5 nsq
#neigh_modify            one 8000

neighbor 0.3 bin

EOF


cat ${GO_PROT1}-out.int_rescale >> in.${GO_PROT1}_lmp


cat >> in.${GO_PROT1}_lmp  <<\EOF

pair_style       hybrid/overlay lj/cut ${bonds_NN} sigmoidal ${cutoff}
pair_coeff * * lj/cut 1 ${sigma}
pair_coeff * * sigmoidal 0 0 0 0
pair_modify shift yes


EOF



cat ${GO_PROT1}-out.paircoeff >> in.${GO_PROT1}_lmp



#insert end commands
cat >> in.${GO_PROT1}_lmp  <<\EOF


velocity      all create 1.0 10

compute cc1 all pair sigmoidal evdwl

thermo_style custom step temp ebond eangle edihed epair pe ke etotal c_cc1
thermo 1000

fix 1 all nve
fix 2 all langevin 1.0 10.0 1 12231 zero yes

dump 1 all atom 1000 unfold.lammpstrj
#dump_modify 1 sort id


thermo          1000
timestep        0.002
run 100000


fix 2 all langevin 10.0 2.0 1 12231 zero yes

thermo          1000
timestep        0.002
run 100000

fix 2 all langevin 2.0 2.0 1 12231 zero yes

thermo          1000
timestep        0.002
run 1000000
EOF
