#!/bin/bash
#BSUB -q gpu_v100
#BSUB -gpu "num=1"
#BSUB -o out.%J
#BSUB -e err.%J

i=2
for j in 10 5 2 1 0.5 0.2 0.1 0.05 0.02 0.01 0
do
                 step="res$(printf "%03d" $i)"
                 nextstep="res$(printf "%03d" $[$i+1])"
                 vmdnumber="$(printf "%03d" $[$i+1])"

                 sed -e "s/MYI/$step/g;s/MYJ/${j}/g;s/MYP/$nextstep/g" ../restraints/restraint_pep_backbone_vmd.tcl > ../restraints/${vmdnumber}_pep_backbone_vmd.tcl
                 echo $i $j

                 vmd -e ../restraints/${vmdnumber}_pep_backbone_vmd.tcl &> ../restraints/${vmdnumber}_pep_backbone_vmd.log
                 wait

                 sed -e "s/MYIN/$step/g" -e "s/MYOUT/$nextstep/g" -e "s/MYRES/${nextstep}_posres.ref/g"  ./res.conf > ./$nextstep.conf

                 namd2 +p5 +devices 0 $nextstep.conf > $nextstep.log
                 [[ $? -ne 0 ]] && break
                 let i=${i}+1
done

i=13
for j in 5 2 1 0.5 0.2 0.1 0.05 0.02 0.01 0
do
                 step="res$(printf "%03d" $i)"
                 nextstep="res$(printf "%03d" $[$i+1])"
                 vmdnumber="$(printf "%03d" $[$i+1])"

                 sed -e "s/MYI/$step/g;s/MYJ/${j}/g;s/MYP/$nextstep/g" ../restraints/restraint_pro_backbone_vmd.tcl > ../restraints/${vmdnumber}_pro_backbone_vmd.tcl
                 echo $i $j

                 vmd -e ../restraints/${vmdnumber}_pro_backbone_vmd.tcl &> ../restraints/${vmdnumber}_pro_backbone_vmd.log
                 wait

                 sed -e "s/MYIN/$step/g" -e "s/MYOUT/$nextstep/g" -e "s/MYRES/${nextstep}_posres.ref/g"  ./res.conf > ./$nextstep.conf

                 namd2 +p5 +devices 0 $nextstep.conf > $nextstep.log
                 [[ $? -ne 0 ]] && break
                 let i=${i}+1
done
