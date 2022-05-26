from optparse import OptionParser
import numpy as np

def removeLowRecombination(in_file_path,recomb_file_path,out_file_path):
    out_file = open(out_file_path,'w')

    rec_intervals=[]
    with open(recomb_file_path) as recomb_file_path:
        while True:
            line = recomb_file_path.readline()
            if not line:
                break
            line_array = line.strip('\n').split('\t')
            rec_intervals.append(list(map(float, line_array[0:2])))
    i=0
    i_max= len(rec_intervals)-1
    end=[x[1] for x in rec_intervals]
    end = np.array(end)
    with open(in_file_path) as in_file:
        while True:
            line = in_file.readline()
            if not line:
                break
            line_array = line.strip('\n').split(',')
            pos=int(line_array[0])
            if any(end<pos): 
                rec_pos=end[end<pos].max()
                i= np.where(end==rec_pos)[0][0]
                if i<i_max-1:
                    if pos>rec_intervals[i][1] and pos < rec_intervals[i+1][0]:
                        out_file.write(line)

def main():
    usage = """%prog  <input> <snp data>"""
    parser = OptionParser(usage)
    parser.add_option("-i", "--inVariantFile", type="string",  default='-',    help="input file path")
    parser.add_option("-r", "--inRecombFile", type="string",  default='-',    help="input file path")
    parser.add_option("-o", "--outFile", type="string",  default='-',    help="output file path")
    options, args= parser.parse_args()
    in_file_path= options.inVariantFile
    recomb_file_path= options.inRecombFile
    out_file_path= options.outFile

    removeLowRecombination(in_file_path,recomb_file_path,out_file_path)


############################## run Main ##############################################################################
if __name__=="__main__":
    main()
