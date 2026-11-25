#!/bin/sh
##
datafile=bto.dat
echo "# u E eta_A_x eta_B_x eta_O_I_x eta_O_II_x eta_O_III_x a b c alpha beta gamma" > $datafile
for s in '1,0,0' '1,1,0' '1,1,1'; do
    ss=`echo $s | sed s/,//g`
    for a in `seq 0.000 0.005 0.121`; do
        amp=`printf "%.3f" $a`
        base="bto.$ss-$amp"
        echo $ss $amp
        sed -e s/0\.045/$amp/ -e s/1,0,0/$s/ valley.py > $base.py
        python $base.py > $base.log
        grep valley: $base.log | sed s/^valley:// >> $datafile
    done
    echo '' >> $datafile
    echo '' >> $datafile
done
