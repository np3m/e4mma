#!/bin/bash
if [[ "$#" -ne 5 ]]; then
    echo "There must be five arguments"
    exit 1
fi
nB=()
Y_e=()
T=()
cs_sq=()

while IFS=, read -r val1 val2 val3
do 
    nB+=($val1)
    Y_e+=($val2)
    T+=($val3)
    cs_temp=$(./eos_nuclei -load $1 -interp-point $val1 $val2 $val3 $5 | awk '/Here:/ {print $2}' | tr -d '\n')
    cs_sq+=($cs_temp)
done < $2

cp $3 $4

for i in ${!cs_sq[@]}; 
do
    acol -read $4 cs2 -value-grid ${nB[i]} ${Y_e[i]} ${T[i]} ${cs_sq[i]} -internal $4
done
