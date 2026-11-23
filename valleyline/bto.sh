#!/bin/sh
##
datafile=bto.dat
echo "# u E eta_A_x eta_B_x eta_O_I_x eta_O_II_x eta_O_III_x a b c alpha beta_ gamma" > $datafile
for a in `seq 0.000 0.005 0.101`; do
    amp=`printf "%.3f" $a`
    base="bto.001.$amp"
    echo $amp
    sed s/0\.045/$amp/ valley.py > $base.py
    python $base.py > $base.log
    grep valley: $base.log | sed s/^valley:// >> $datafile
done
