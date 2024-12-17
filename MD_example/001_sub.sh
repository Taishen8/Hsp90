#!/bin/bash
#BSUB -q gpu_v100
#BSUB -gpu "num=1"
#BSUB -o out.%J
#BSUB -e err.%J

cd restraints
vmd -e 001_minimize_vmd.tcl
cd ..
wait
cd res
namd2 +p5 +devices 0 res001.conf > res001.log
cd ..
wait
cd restraints
vmd -e 002_nvt_vmd.tcl
cd ..
wait
cd res
namd2 +p5 +devices 0 res002.conf > res002.log 
wait
nohup bash 001_sub_res.sh &
cd ..
wait
cd equ
nohup bash 001_sub_equ.sh &
