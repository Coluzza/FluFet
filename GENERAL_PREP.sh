#!/bin/bash
#it reads the pdb identifies each chain but still prints out a single object for the simulation. it will also give you the length of each chain so that you can easily write the ignore.dat file
PDB="$1"
PARAM="$2"
pka="$3"
pH="$4"

chi=$(basename $PDB | awk -F ".pdb" '{print $1}')

CHAINS=$(awk '{if((substr($0,1,6)=="ATOM  ")&&(substr($0,14,2)=="CA")) print substr($0,22,1)}' "$1" | sort -u)
NUM_CHAINS=$(echo "$CHAINS" | wc -w)
echo "Chains: $CHAINS"
echo $NUM_CHAINS
top=$(grep -n SUMMARY  $pka | awk -F ":" '{print $1}')
bottom=$(grep -n "Free energy"  $pka | awk -F ":" '{print $1}') 

awk -v t=$top -v b=$bottom 'NR>t+1&&NR<b-2&&$1!="N+"&&$1!="C-"{print $2,$1,$4}' $pka | sort -k 1n,1 >"$chi"_summary.pka
pka_mask.awk "pH=$pH" "$chi"_summary.pka > temp
mv temp "$chi"_summary.pka
if [ -f grofiles.list ]; then
    rm grofiles.list
fi


length=()
j=0
for chain in $CHAINS; do
    echo "Processing chain: $chain"
    awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,22,1)=="'$chain'")&&(substr($0,18,3)!="HOH")&&((substr($0,14,1)=="S")||(substr($0,14,2)=="CB")||(substr($0,14,2)=="CA")||(substr($0,14,2)=="C ")||(substr($0,14,2)=="N ")||(substr($0,14,2)=="O ")))||(substr($0,1,6)=="ENDMDL")) print $0}' "$1" | awk '{if(($1=="ENDMDL")) exit; else print $0}' > temp
    end=$(awk '{if(substr($0,14,2)=="O ") print NR}' temp | tail -n 1)
    start=$(awk '{if(substr($0,14,2)=="N ") print NR}' temp | head -n 1)
    echo "$start $end"
    awk '{if((NR>='$start')&&(NR<='$end')) print substr($0,14,2),substr($0,18,3),substr($0,31,8),substr($0,39,8),substr($0,47,8),substr($0,17,1)}' temp > "temp_CMAPHOH.pdb"

    length+=($(grep CA "temp_CMAPHOH.pdb" | wc -l))

    rm temp temp_CMAPHOH.pdb
    echo "$chain ${length[$j]}"

    ((j++))
done

if (( $NUM_CHAINS > 1 )); then
    nignore=$(echo $NUM_CHAINS | awk '{print $1-1}')
    echo "$nignore" > "ignore_${chi}.dat"
    j=0
    cumulative_length=0
    while (( $j < $nignore )); do
        cumulative_length=$((cumulative_length + ${length[$j]}))
        echo "$cumulative_length" >> "ignore_${chi}.dat"
        ((j++))
    done
fi

#awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,22,1)=="'$chain'")&&(substr($0,18,3)!="HOH")&&((substr($0,14,1)=="S")||(substr($0,14,2)=="CB")||(substr($0,14,2)=="CA")||(substr($0,14,2)=="C ")||(substr($0,14,2)=="N ")||(substr($0,14,2)=="O ")))||(substr($0,1,6)=="ENDMDL")) print $0}' "$1" | awk '{if(($1=="ENDMDL")) exit; else print $0}' > temp
awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,18,3)!="HOH")&&((substr($0,14,1)=="S")||(substr($0,14,2)=="CB")||(substr($0,14,2)=="CA")||(substr($0,14,2)=="C ")||(substr($0,14,2)=="N ")||(substr($0,14,2)=="O ")))||(substr($0,1,6)=="ENDMDL")) print $0}' "$1" | awk '{if(($1=="ENDMDL")) exit; else print $0}' > temp

end=$(awk '{if(substr($0,14,2)=="O ") print NR}' temp | tail -n 1)
start=$(awk '{if(substr($0,14,2)=="N ") print NR}' temp | head -n 1)
echo "$start $end"

awk '{if((NR>='$start')&&(NR<='$end')) print substr($0,14,2),substr($0,18,3),substr($0,31,8),substr($0,39,8),substr($0,47,8),substr($0,17,1)}' temp > "${chi}_CMAPHOH.pdb"
awk '{if(((substr($0,1,6)=="ATOM  ")&&(substr($0,18,3)!="HOH")&&((substr($0,14,2)=="CA")))||(substr($0,1,6)=="ENDMDL")) print substr($0,23,4),substr($0,18,3)}' temp | awk '{print $1,$2}'  > "$chi"_MASK.pdb

rm temp




seq_mask.awk "$chi"_MASK.pdb > temp
cat "$chi"_summary.pka temp | sort -k 1n,1 > "$chi"-out.pka
pka_to_mask.awk "$chi"-out.pka > "$chi"-charge.mask 


echo "$chi"
sed  "s/Protein_PDB.*#/Protein_PDB = ${chi}_CMAPHOH.pdb #/g" $PARAM > temp1
sed  "s/Ignore_Bond_List.*#/Ignore_Bond_List = ignore_${chi}.dat #/g" temp1 > temp2
sed  "s/Prefix.*#/Prefix = ${chi} #/g" temp2 > temp1
sed  "s/Protein_Mask_Protonation.*#/Protein_Mask_Protonation = "$chi"-charge.mask #/g" temp1> param.dat

GraftFlow_Brush_general.x
