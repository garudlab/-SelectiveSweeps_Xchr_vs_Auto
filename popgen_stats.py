from optparse import OptionParser
import numpy as np
import math
import linecache


def get_popgen_stats(in_file_path,out_file_path,out_file2_path,N):
    in_file =  open(in_file_path,'r')
    out_file = open(out_file_path,'w')
    out_file2 = open(out_file2_path,'a')
    out_file.write('position'+','+'S'+","+'pi'+","+'watersons'+"\n")
    #out_file.write('position'+','+'S'+","+'pi'+","+"\n")
    binLength=round(N, -4)

    #binLength=int(math.ceil(N / 1e4)) * 1e4
    step=1e4
    bins=np.arange(0,binLength+step,step)

    line=in_file.readline()
    snp_list = line.strip('\n').split(',')
    pos= int(snp_list[0])
    snps10kb=[]
    heterozigosity_lst=[]
    pi=0
    pi_Full=0
    for i in range(1,len(bins)):
        if pos<N:
            while pos < bins[i]:
                snps10kb.append(snp_list)
                line=in_file.readline()
                snp_list = line.strip('\n').split(',')
                if snp_list == ['']:
                    break;
                else:
                    pos= int(snp_list[0])
                    #H=heterozigosity(snp_list)
                    pi=pi+heterozigosity2(snp_list)
                    pi_Full=pi_Full+heterozigosity2(snp_list)
                    #heterozigosity_lst.append(H)
            out_file.write(str(bins[i])+','+str(len(snps10kb)/1e4)+","+str(pi/1e4)+","+str(len(snps10kb)/(5.550497*1e4))+"\n")
            #out_file.write(str(bins[i])+','+str(len(snps10kb)/1e4)+","+str(sum(heterozigosity_lst)/1e4)+","+str(len(snps10kb)/(5.550497*1e4))+"\n")
            #calculate Sper10kb
            #S_per_10kb(i,out_file,bins,snps10kb)
            #calculate pi per 10Kb
            #pi_per_10kb(i,out_file,bins,snps10kb)
            snps10kb=[]
            #heterozigosity_lst=[]
            pi=0
    out_file2.write(str(pi_Full)+"\n")
    
    in_file.close()
    out_file.close()

def heterozigosity(snp_list):
    snps=snp_list[1:]
    #get snp counts
    snp_dict = {key:snps.count(key) for key in snps}
    #get frequencies
    snp_dict =  {key: value / total for total in (sum(snp_dict.values()),) for key, value in snp_dict.items()}
    #remove N's
    snp_dict = { key: value for key, value in snp_dict.items() if key != "N" }
    #square values of dictionary
    snp_dict= { key: value**2 for key, value in snp_dict.items()}
    #get sum of suares
    G=sum(snp_dict.values())
    return 1-G

def heterozigosity2(snp_list):
    snps=snp_list[1:]
    #get snp counts
    unique_bases = set(snps)
    if ('N' in unique_bases): unique_bases.remove('N')
    if (len(unique_bases)==2):
        snp_dict = {key:snps.count(key) for key in snps}
        #remove N's
        snp_dict = { key: value for key, value in snp_dict.items() if key != "N" }
        sample_size = sum(snp_dict.values())
        #get frequencies
        snp_dict =  {key: value / total for total in (sum(snp_dict.values()),) for key, value in snp_dict.items()}
        p=max(snp_dict.values())
    else:
        p=0
        sample_size = 2
    return 2*p*(1-p)*sample_size/(sample_size-1)

def S_per_10kb(i,out_file,bins,snps10kb):
    out_file.write(str(bins[i])+','+str(len(snps10kb)/1e4)+"\n")

def pi_per_10kb(i,out_file,bins,snps10kb):
    if len(snps10kb)>0: print(snps10kb[0])

def get_LD(in_file_path,out_file_path,N):

    in_file =  open(in_file_path,'r')
    out_file = open(out_file_path,'w')

    numberLines =  len(in_file.readlines())
    jump=50
    window=1e4

    snp_list = linecache.getline(in_file_path, 1).strip().split(',')

    posStart=int(snp_list[0])
    pos=posStart
    lastSNP = int(linecache.getline(in_file_path, numberLines).split(',')[0]) - window

    majorAllele_lst = []
    snps10kb=[]
    i=1
    iter_indx=i+1
    posStart_new = int(linecache.getline(in_file_path, iter_indx).strip().split(',')[0])
    
    #print(lastSNP)
    while pos <= lastSNP:
        print(posStart)
        snp_list = linecache.getline(in_file_path, i).strip().split(',')
        pos=int(snp_list[0])
        if pos - posStart<=window:
            if pos - posStart<=50 and pos-posStart >0: #find start position when iterating by 50 bp
                posStart_new = pos
                iter_indx = i

            majorAl = get_major_allele(snp_list)
            if majorAl != 'Null':
                majorAllele_lst.append(majorAl)
                snps10kb.append(snp_list)
            i+=1

        elif len(majorAllele_lst)>=4 :
            [LD_pairwise_array, distance_array]=get_pairwiseLD_Rsquared(snps10kb,majorAllele_lst)
            for r in range(0, len(LD_pairwise_array)):
                out_file.write(str(LD_pairwise_array[r]) + '\t' + str(distance_array[r]) + '\n')
           
            posStart=posStart_new
            i=iter_indx
            iter_indx +=1
            posStart_new = int(linecache.getline(in_file_path, iter_indx).strip().split(',')[0])

            majorAllele_lst = []
            snps10kb=[]

        else:
            posStart=posStart_new
            i=iter_indx
            iter_indx +=1
            posStart_new = int(linecache.getline(in_file_path, iter_indx).strip().split(',')[0])

            

def get_major_allele(snp_list):
    snps=snp_list[1:]
    snp_dict = {key:snps.count(key) for key in snps}
    #remove N's
    snp_dict = { key: value for key, value in snp_dict.items() if key != "N" }

    majorAllele = 'Null'
    if (len(snp_dict)==2):
        snp_dict= { key: value/float(numStrains) for key, value in snp_dict.items()}
        if max(snp_dict.values())<0.95 and max(snp_dict.values())>0.05:
            majorAllele = max(snp_dict, key=snp_dict.get)

    return majorAllele

def get_pairwiseLD_Rsquared(snps10kb,majorAllele_lst):
    LD_pairwise_array =[]
    distance_array = []
    for i in range(1,len(snps10kb)):
        A= majorAllele_lst[0]
        B = majorAllele_lst[i]
        # Iterate again through snps10kb and calculate P(AB)
        SNPAB = 0
        SNPAB_denom = 0
        SNPA = 0
        SNPA_denom = 0
        SNPB = 0
        SNPB_denom = 0
        #if len(snps10kb)==1 : print(snps10kb)
        for j in range(1,numStrains+1):
            if snps10kb[0][j] == A and snps10kb[i][j]==B:
                SNPAB += 1
                SNPAB_denom += 1
                SNPA += 1
                SNPA_denom += 1
                SNPB += 1
                SNPB_denom += 1
            elif snps10kb[0][j] !='N':
                if snps10kb[0][j] == A:
                    SNPA +=1
                if snps10kb[i][j] == B:
                    SNPB +=1
                SNPAB_denom += 1
                SNPA_denom += 1
                SNPB_denom += 1
        if SNPA_denom>0:
            P_a = float(SNPA)/SNPA_denom
        else:
            P_A=0
        if SNPB_denom >0:
            P_b =float(SNPB)/SNPB_denom
        else:
            P_b=0
        if SNPA != SNPA_denom and SNPB != SNPB_denom:
            P_ab = float(SNPAB)/(SNPAB_denom)

            if P_a!=0 and P_b!=0 and P_a!=1 and P_b!=1: # The reason for this check is that allele frequency can go to zero when checking for Ns in the other SNP.
                LD_pairwise_array.append(((P_ab - P_a*P_b)**2)/(P_a*(1-P_a)*P_b*(1-P_b)))
                distance_array.append(abs(int(snps10kb[i][0])-int(snps10kb[0][0])))

    return [LD_pairwise_array, distance_array]

def main():
    usage = """%prog  <input> <snp data>"""
    parser = OptionParser(usage)
    parser.add_option("-i", "--inFile", type="string",  default='-',    help="input file path for output")
    parser.add_option("-o", "--outFile", type="string",  default='-',    help="output file path")
    parser.add_option("-a", "--outFile2", type="string",  default='-',    help="output file2 path")
    parser.add_option("-N", "--nLines", type="int",  default=1,    help="last position in file")
    options, args= parser.parse_args()

    in_file_path= options.inFile
    out_file_path= options.outFile
    out_file2_path= options.outFile2
    N = options.nLines

    global numStrains
    numStrains = int(args[0])
    #get_popgen_stats(in_file_path,out_file_path,out_file2_path,N)
    get_LD(in_file_path,out_file_path,N)


############################## run Main ##############################################################################
if __name__=="__main__":
    main()
