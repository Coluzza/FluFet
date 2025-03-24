#!/usr/bin/awk -f


$2=="ALA"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="CYS"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="ASP"{if(pH>$3) print $1,$2,$3,0,1,-1; else print $1,$2,$3,0,1,0}
$2=="GLU"{if(pH>$3) print $1,$2,$3,0,1,-1; else print $1,$2,$3,0,1,0}
$2=="PHE"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="GLY"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="HIS"{if(pH>$3) print $1,$2,$3,1,1,0; else print $1,$2,$3,1,1,1}
$2=="ILE"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="LYS"{if(pH>$3) print $1,$2,$3,1,0,0; else print $1,$2,$3,1,0,1}
$2=="LEU"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="MET"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="ASN"{if(pH>$3) print $1,$2,$3,1,1,0; else print $1,$2,$3,1,1,0}
$2=="PRO"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="GLN"{if(pH>$3) print $1,$2,$3,1,1,0; else print $1,$2,$3,1,1,0}
$2=="ARG"{if(pH>$3) print $1,$2,$3,1,0,0; else print $1,$2,$3,1,0,1}
$2=="SER"{if(pH>$3) print $1,$2,$3,1,1,0; else print $1,$2,$3,1,1,0}
$2=="THR"{if(pH>$3) print $1,$2,$3,1,1,0; else print $1,$2,$3,1,1,0}
$2=="VAL"{if(pH>$3) print $1,$2,$3,0,0,0; else print $1,$2,$3,0,0,0}
$2=="TRP"{if(pH>$3) print $1,$2,$3,1,0,0; else print $1,$2,$3,1,0,0}
$2=="TYR"{if(pH>$3) print $1,$2,$3,1,0,0; else print $1,$2,$3,1,0,0}

