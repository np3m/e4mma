#!/bin/bash
# $1 csv of superluminal points with fixed speed of sound.
# $2 st.o2
# $3 Where fixed st.o2 file will be saved
if [[ "$#" -ne 3 ]]; then
    echo "There must be three arguments"
    exit 1
fi

cp $2 $3

while IFS=, read -r val1 val2 val3 val4
do
    acol -read $3 cs2 -value-grid $val1 $val2 $val3 $val4 -internal $3
done < $1
