#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file.
file=$1
rho_in=$2
PF_on=$3
Ne=$4
Q=$5
sd=$6
sb=$7
chr=$8
burn_in=$9
H=${10}
ID=${11}
sexRatio=0.5
mu=1e-9

#outFile=outFiles/ConstantNe${Ne}_theta_${adaptive_theta}_selection_${selection}_Q${Q}_Chr${chr}
#outFile=outFiles/ConstantNe_theta_selection_${selection}_Chr${chr}

for i in `seq 1 100`; do #1100

        echo $id
        #run SLiM
	slim -d "file_path_data='${file}_${ID}'" -d "file_path='${file}_${ID}.trees'" -d R=${rho_in} -d PF_on=${PF_on} -d N=${Ne} -d H=${H} -d run=${ID}$i \
        -d Q=${Q} -d sd=${sd} -d sb=${sb} -d "ChrType='${chr}'" -d burn_in=${burn_in} \
        -d sexRatio=${sexRatio} sweeps_SGV_clusterHaps.slim

    #add mutations and recapitate
        python3.7 mutate_and_recapitate_SGVhaps.py ${file}_${ID}.trees ${file}_numClusterHaps_$ID $rho_in $Ne $mu $Q ${PF_on} $H $chr

	rm ${file}_${ID}.trees
done
