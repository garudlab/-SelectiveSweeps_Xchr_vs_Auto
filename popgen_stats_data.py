from optparse import OptionParser
import numpy as np
import math
import linecache

def get_popgen_stats(in_file,out_file,w):

    out_file = open(out_file,'a')
    #out_file2 = open(out_file2_path,'w')
    out_file.write('position'+'\t'+'S'+"\t"+'pi'+"\t"+'watersons'+"\n")

    countLines =  open(in_file)
    numberLines =  len(countLines.readlines())
    last_pos = int(linecache.getline(in_file,numberLines).split(',')[0].strip('\n'))
    last_pos_bins= round(last_pos,-4)

    #binLength=int(math.ceil(N / 1e4)) * 1e4
    bins=np.arange(0,last_pos_bins+w*1000,w*1000)
    snps10kb=[]
    heterozigosity_lst=[]
    
    line=linecache.getline(in_file,1).strip('\n').split(',')
    pos = int(line[0])
    snp_list=line[1:]
    first_bin=bins[bins>=pos][0] # find first bin
    j=np.absolute(bins-first_bin).argmin() #find index of first bin in bins array
    pi=0
    pi_Full=0
    i=1
    for window in bins[j:]: #calc pi in non-overlapping 10kb windows
        #print(window)
        if pos<=last_pos:
            while pos < window:
                snps10kb.append(snp_list)
                pi=pi+heterozigosity2(snp_list)
                pi_Full=pi_Full+heterozigosity2(snp_list)
                line=linecache.getline(in_file,i).strip('\n').split(',')
                if line == ['']:
                    break;
                else:
                    pos = int(line[0])
                    snp_list=line[1:]
                i+=1
            out_file.write(str(window)+'\t'+str(len(snps10kb)/1e4)+"\t"+str(pi/1e4)+"\t"+str(len(snps10kb)/(5.550497*1e4))+"\n")
            snps10kb=[]
            pi=0  

    #out_file2.write(str(pi_Full)+"\n")
    out_file.close()

def heterozigosity2(snp_list):
    #get snp counts
    unique_bases = set(snp_list)
    if ('N' in unique_bases): unique_bases.remove('N')
    if (len(unique_bases)==2):
        snp_dict = {key:snp_list.count(key) for key in snp_list}
        snp_dict_all = {key: value / total for total in (sum(snp_dict.values()),) for key, value in snp_dict.items()}
        snp_dict_all =  { key: value for key, value in snp_dict_all.items() if key != "N" } #remove N
        #remove N's
        snp_dict = { key: value for key, value in snp_dict.items() if key != "N" }
        #get frequencies
        snp_dict =  {key: value / total for total in (sum(snp_dict.values()),) for key, value in snp_dict.items()}
        #p=max(snp_dict.values())
        p=max(snp_dict_all.values())
        q=min(snp_dict_all.values())
    else:
        p=0
    return 2*p*q*numStrains/(float(numStrains-1))

def main():
    usage = """%prog  <input> <snp data>"""
    parser = OptionParser(usage)
    parser.add_option("-i", "--inFile", type="string",  default='-',    help="input file path for output")
    parser.add_option("-o", "--outFile", type="string",  default='-',    help="output file path")
    #parser.add_option("-a", "--outFile2", type="string",  default='-',    help="output file2 path")
    parser.add_option("-w", "--window", type="int",  default=10,    help="window length in kb")
    options, args= parser.parse_args()

    in_file= options.inFile
    out_file= options.outFile
    #out_file2_path= options.outFile2
    w = options.window

    global numStrains
    numStrains = int(args[0])
    get_popgen_stats(in_file,out_file,w)
    #get_LD(in_file_path,out_file_path,N)

############################## run Main ##############################################################################
if __name__=="__main__":
    main()