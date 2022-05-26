#!/bin/bash
# This script takes in an MS file  with 1000 simulations, and outputs a concatenated processed file. 
file=$1
rho_in=$2
adaptive_theta=$3
selection=$4
locusLength=$5
NAm=$6
EUr=$7


outFile=outFiles/Admixture_fixedPopSize_NAm${NAm}_EUr${EUr}_theta_${adaptive_theta}_selection_${selection}_MS

for i in `seq 1 2500`; do #`seq 1 1000`
    #echo $i
    python admixture_parameters_fixedPopSize_MS.py $rho_in $adaptive_theta $selection $locusLength $NAm $EUr > ${file}_var

    command=`cat ${file}_var | head -1` 
    Nac=`cat ${file}_var | head -3 | tail -1` 
    s=`cat ${file}_var | head -4 | tail -1` 
    age=`cat ${file}_var | head -5 | tail -1` 


    eval $command

    PF=`cat $file | grep -B 1 segsites | head -1 | cut -f7`
    segsites=`cat $file | grep segsites | cut -f2 -d' '`
    echo 'PF is' $PF
    echo 'segsites is' $segsites
    echo 'Nac is' $Nac
    echo 's is' $s
    echo 'age is' $age

# proceed with analyzing the sweep                                   

# cut the file

lineNo=`cat $file | grep -n segsites | cut -f1 -d ':'`

(( lineNo = lineNo + 101)) #146
cat $file | head -${lineNo} | tail -101 > ${file}_cut #146


numLines=`cat ${file}_cut | egrep '(2|3|4|5|6|7|8|9|a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z)' | wc -l`

if [ $numLines == 1 ]
    then
        HS=1
    else
        HS=0

fi


#convert to MSMS format

if [ "$segsites" = 0 ]
    then
    python convertMS_noSegSites.py 100 ${file}_MS
    else
    python convertMS.py ${file}_cut ${file}_MS
    python removeSingletons.py ${file}_MS ${file}_noSingletons
fi


# cluster 


if [ "$locusLength" -eq "100000" ]; then

     if [ "$selection" = True ]
     then
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_401 -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age -m $segsites 
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_265 -w 265 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age -m $segsites 
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age -m $segsites -z 265

     else
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_401 -w 401 -j 50 -d 0 -s 0.5
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_265 -w 265 -j 50 -d 0 -s 0.5
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -z 265

     fi

     cat ${file}_cluster_snps_401 >> ${outFile}_snps_401.txt
     cat ${file}_cluster_snps_265 >> ${outFile}_snps_265.txt
     cat ${file}_cluster_snps >> ${outFile}_snps.txt

fi





if [ "$locusLength" -eq "600000" ]; then
     if [ "$selection" = True ]
     then
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_401 -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age -m $segsites 
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_265 -w 265 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age -m $segsites 
 	    python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -p $PF -n $Nac -e $s -a $age -m $segsites -z 265

     else
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_401 -w 401 -j 50 -d 0 -s 0.5
        python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps_265 -w 265 -j 50 -d 0 -s 0.5
 	    python H12_H2H1_simulations.py ${file}_MS 100 -o ${file}_cluster_snps -w 401 -j 50 -d 0 -s 0.5 -z 265

     fi

     cat ${file}_cluster_snps_401 >> ${outFile}_snps_401.txt
     cat ${file}_cluster_snps_265 >> ${outFile}_snps_265.txt
     cat ${file}_cluster_snps >> ${outFile}_snps.txt

 fi



rm ${file}
rm ${file}_var
rm ${file}_cut
rm ${file}_MS
rm ${file}_cluster_snps*
rm ${file}_noSingletons
##rm ${file}_cluster_bps

done




