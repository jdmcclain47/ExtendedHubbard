#!/bin/bash

prog="a.out"
outfile="TEMP.OUT"

for i in gamma_6_6_6 gamma_7_7_7
#gamma_2_2_2 gamma_3_3_3 gamma_4_4_4 gamma_5_5_5 #
do
	mkdir $i;
	arr=(${i//_/ });
        echo "${arr[2]}"
for j in WS #EWALD EWALDMSE
do
	cd $i; mkdir $j; cd $j; echo "${arr[1]} ${arr[2]} ${arr[3]}" > KPOINT;
        cd ../../;
	cp INCAR $i/$j;
	cp PARAM $i/$j;
	cp POSCAR $i/$j;
        echo "XCMETHOD $j" >> $i/$j/PARAM; 
	cd $i; cd $j;
	#nohup "./../../$prog" > $outfile &
        cd ../../;
done
done
