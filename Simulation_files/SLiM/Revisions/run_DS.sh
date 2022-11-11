#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 
file=$1
rho_in=$2
adaptive_theta=$3
Ne=$4
Q=$5
sd=$6
sb=$7
chr=$8
burn_in=$9
H=${10}
sexRatio=0.5

for i in `seq 1 20`; do #1100

	#runSlim
	slim -d "file_path_data='${file}'" -d R=${rho_in} -d THETA_A=${adaptive_theta} -d N=${Ne} -d H=${H} \
    -d Q=${Q} -d sd=${sd} -d sb=${sb} -d "ChrType='${chr}'" -d burn_in=${burn_in} -d sexRatio=${sexRatio} DominanceShifts_hbFix.slim

done
