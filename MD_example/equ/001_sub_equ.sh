#!/bin/bash
#BSUB -q gpu_v100
#BSUB -gpu "num=1"
#BSUB -o out.%J
#BSUB -e err.%J

#cp ../res/res023.coor ./
#cp ../res/res023.vel ./
#cp ../res/res023.xsc ./
#mv res023.coor equ000.coor
#mv res023.vel equ000.vel
#mv res023.xsc equ000.xsc

for i in {361..499}
do
        step="equ$(printf "%03d" $i)"
        nextstep="equ$(printf "%03d" $[$i+1])"
        sed -e "s/MYIN/$step/g" -e "s/MYOUT/$nextstep/g"  ./equ.conf > ./$nextstep.conf
        echo $i
        echo $j
        namd2 +p44 +devices 0 $nextstep.conf > $nextstep.log
        [[ $? -ne 0 ]] && break
done
