#!/bin/bash
#set -x

FILEIN="$1"
FILEOUT="$2"

echo $FILEIN
echo $FILEOUT

#count output lines
wc -l $FILEIN > tempLine.$FILEIN
linecount=$(awk '{print $1}' tempLine.$FILEIN)
echo $linecount
rm tempLine.$FILEIN

#remove 3 lines of header
headremove=$(bc <<< "scale=2;$linecount-3")
echo $headremove
tail -n $headremove $FILEIN > tempHead.$FILEIN

#count again
wc -l tempHead.$FILEIN > tempLine.$FILEIN
linecount=$(awk '{print $1}' tempLine.$FILEIN)
echo $linecount
rm tempLine.$FILEIN

cutnumber1=$(bc <<< "scale=0;$linecount/2-1")
cutnumber2=$(bc <<< "scale=0;$linecount/2")


tail -n $cutnumber1 tempHead.$FILEIN > $FILEOUT.1
head -n $cutnumber2 tempHead.$FILEIN > temptempSplit.$FILEIN.2
tail -n $cutnumber1 temptempSplit.$FILEIN.2 > $FILEOUT.2
rm temptempSplit.$FILEIN.2
rm tempHead.$FILEIN




#rm tempHead.$FILEIN



#skip=$(bc <<< "scale=2;$N+9")

