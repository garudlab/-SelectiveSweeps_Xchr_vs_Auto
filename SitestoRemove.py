

######################read files###############################################
missing_perSite_path =  open("missingDat_perSite.txt",'r')
out_file = open("toRemove.txt",'w')

missing_perSite = missing_perSite_path.read().splitlines()
###############################################################################

missing_perSite = list(map(int, missing_perSite))
missing_perSite = [i / 100 for i in missing_perSite]
removeLines=[ n for n,i in enumerate(missing_perSite) if i>0.5 ]

[out_file.write(str(i+1)+'\n') for i in removeLines]
###############################################################################
