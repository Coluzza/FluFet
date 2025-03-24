#!/bin/bash

pH=$3
chi=$(echo $1 | awk -F ".pdb" '{print $1}')
CHAINS=$(awk '{if((substr($0,1,6)=="ATOM  ")&&(substr($0,14,2)=="CA")) print substr($0,22,1)}' $1 | sort -u)


top=$(grep -n SUMMARY  $2 | awk -F ":" '{print $1}')
bottom=$(grep -n "Free energy"  $2 | awk -F ":" '{print $1}') 

awk -v t=$top -v b=$bottom 'NR>t+1&&NR<b-2&&$1!="N+"&&$1!="C-"{print $2,$1,$4}' $2 | sort -k 1n,1 >"$chi"_summary.pka
pka_mask.awk "pH=$pH" "$chi"_summary.pka > temp
mv temp "$chi"_summary.pka
if [ -z $CHAINS ]
then

awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,18,3)!="HOH")&&((substr($0,14,1)=="S")||(substr($0,14,2)=="CB")||(substr($0,14,2)=="CA")||(substr($0,14,2)=="C ")||(substr($0,14,2)=="N ")||(substr($0,14,2)=="O ")))||(substr($0,1,6)=="ENDMDL")) print $0}' $1 | awk '{if($1=="ENDMDL") exit; else print $0}' > temp
skip=$(awk '{if(substr($0,14,2)=="O ") print $2}' temp | tail -n 1)

awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,18,3)!="HOH")&&((substr($0,14,2)=="CA")))||(substr($0,1,6)=="ENDMDL")) print substr($0,18,3)}' temp > "$chi"_CMAPHOH.pdb
rm temp
echo $chi

seq_mask.awk "$chi"_CMAPHOH.pdb

else
for  chain in $CHAINS
do
awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,22,1)=="'$chain'")&&(substr($0,18,3)!="HOH")&&((substr($0,14,1)=="S")||(substr($0,14,2)=="CB")||(substr($0,14,2)=="CA")||(substr($0,14,2)=="C ")||(substr($0,14,2)=="N ")||(substr($0,14,2)=="O ")))||(substr($0,1,6)=="ENDMDL")) print $0}' $1 | awk '{if($1=="ENDMDL") exit; else print $0}' > temp
skip=$(awk '{if(substr($0,14,2)=="O ") print $2}' temp | tail -n 1)

awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,18,3)!="HOH")&&((substr($0,14,2)=="CA")))||(substr($0,1,6)=="ENDMDL")) print substr($0,23,4),substr($0,18,3)}' temp | awk '{print $1,$2}'  > "$chi"_"$chain"_CMAPHOH.pdb

seq_mask.awk "$chi"_"$chain"_CMAPHOH.pdb > temp
cat "$chi"_summary.pka temp | sort -k 1n,1 > "$chi"_"$chain".pka
pka_to_mask.awk "$chi"_"$chain".pka > "$chi"_"$chain".mask 

#while read c
#do
#pka=$(echo $c | awk '{print $3}')
#idx=$(echo $c | awk '{print $1}')
#awk -v pka=$pka -v idx=$idx '{if($1==idx) print $0,pka}' temp 
#done < "$chi"_summary.pka

rm temp
echo "$chi"_$chain


done

fi
