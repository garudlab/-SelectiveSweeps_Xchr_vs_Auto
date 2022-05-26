#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 
file=$1
rho_in=$2
adaptive_theta=$3
selection=$4
locusLength=$5
Ne=$6
mu=1e-9
id=$7
chr=$8
Q=$9


outFile=outFiles/ConstantNe${Ne}_theta_${adaptive_theta}_selection_${selection}_Q${Q}_Chr${chr}
#outFile=outFiles/ConstantNe_theta_selection_${selection}_Chr${chr}

for i in `seq 1 300`; do #1100 

	echo $id
	#get slim command
	python ConstantNe_parameters.py $rho_in $adaptive_theta $selection $Ne $Q $chr slim $id $i> ${file}_var

	command=`cat ${file}_var | head -1` 

	echo $command

    eval $command

    PF=`echo $command | grep -o "PF=.*" | tr ' ' '\n' | head -1 |cut -d"=" -f2-`
    s=`echo $command | grep -o "sb=.*" | tr ' ' '\n' | head -1| cut -d"=" -f2- `
    age=`echo $command | grep -o "AgeSweep=.*" | tr ' ' '\n' | head -1| cut -d"=" -f2-`
    echo 'PF is' $PF
    echo 's is' $s
    echo 'age is' $age

    #add mutations and recapitate
    python3.7 mutate_and_recapitate.py tmp_intermediate_files/ConstantNe_${id}${i}.trees ${file}_vcf $rho_in $Ne $mu $Q ${file}_haps
    #i can also try oputputting haplotypes directly

    #parse VCF for H12
    #python Parse_SlimVCF.py -i ${file}_vcf -o ${file}_vcf_Parsed
    #remove singletons
    #python removeSingletons.py ${file}_vcf_Parsed ${file}_vcf_noSingletons
    python removeSingletons.py ${file}_haps ${file}_haps_noSingletons

    

    segsites=`cat ${file}_haps |wc -l`
    #run H12

    echo 'num snps is' $segsites

    if [ "$selection" = True ]
     then
        python H12_H2H1_simulations.py ${file}_haps 100 -o ${file}_cluster_snps_401 -w 401 -j 50 -d 0 -s 50000 -p $PF -n $Ne -e $s -a $age -m $segsites -q $Q
        python H12_H2H1_simulations.py ${file}_haps 100 -o ${file}_cluster_snps_265 -w 265 -j 50 -d 0 -s 50000 -p $PF -n $Ne -e $s -a $age -m $segsites -q $Q
 	python H12_H2H1_simulations.py ${file}_haps 100 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 50000 -p $PF -n $Ne -e $s -a $age -m $segsites -z 265 -q $Q

     else
        python H12_H2H1_simulations.py ${file}_haps 100 -o ${file}_cluster_snps_401 -w 401 -j 50 -d 0 -s 50000 -m $segsites  -x $adaptive_theta -q $Q
        python H12_H2H1_simulations.py ${file}_haps 100 -o ${file}_cluster_snps_265 -w 265 -j 50 -d 0 -s 50000 -m $segsites  -x $adaptive_theta -q $Q
 	python H12_H2H1_simulations.py ${file}_haps 100 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 50000 -z 265 -m $segsites  -x $adaptive_theta -q $Q

     fi

     cat ${file}_cluster_snps_401 >> ${outFile}_snps_401.txt
     cat ${file}_cluster_snps_265 >> ${outFile}_snps_265.txt
     cat ${file}_cluster_snps >> ${outFile}_snps.txt


rm ${file}_var
#rm ${file}_vcf
#rm ${file}_vcf_Parsed
rm ${file}_haps
rm ${file}_cluster_snps*
#rm ${file}_vcf_noSingletons
rm ${file}_haps_noSingletons
rm tmp_intermediate_files/ConstantNe_${id}${i}.trees
done
