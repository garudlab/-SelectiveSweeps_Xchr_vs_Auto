#! /usr/bin/python
# bedTileElements - tile each record in input.bed with potentially shorter elements which are written to output.bed 
import sys
from optparse import OptionParser
import copy
import os

####################################################################################
# path psyco
####################################################################################
sys.path.append("/home/jsp/prog/utillities/py_modules")
# to speed things up
#import psyco
#psyco.full()
#import bed

###############################


def reformatForJoin(inFile, outFile):

    for line in inFile:
        line=line.strip('\n')
        genotypes=line.split(',')
        nucleotides=[0,0,0,0]
        for i in range(1, len(genotypes)):
            if genotypes[i] == 'A':
                nucleotides[0] +=1
            if genotypes[i] == 'T':
                nucleotides[1] +=1
            if genotypes[i] == 'G':
                nucleotides[2] +=1
            if genotypes[i] == 'C':
                nucleotides[3] +=1
        # check whether or not this SNP has two or more alleles
        counter=0
        for y in range(0, len(nucleotides)):
            if nucleotides[y]>0:
                counter +=1
        if  counter >= 2:
            outFile.write(','.join(genotypes[0:(len(genotypes)-1)]) + '\n')



#######################
def mkOptionParser():
    """ Defines options and returns parser """
    
    usage = """%prog  <input> <output>
    %prog removes the invariant sites in the DGRP file after removeing the bad haplotypes and replacing missing data with Ns"""

    parser = OptionParser(usage)
   

    return parser



#####################


def main():
    """ see usage in mkOptionParser. """
    parser = mkOptionParser()
    options, args= parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments")


    inFN         = args[0]
    outFN        = args[1]
    
    if inFN == '-':
        inFile = sys.stdin
    else:
        inFile      = open(inFN, 'r')
   
    if outFN == '-':
        outFile = sys.stdout
    else:
        outFile      = open(outFN, 'w')



    reformatForJoin(inFile, outFile)


    

#run main
if __name__ == '__main__':
    main()

