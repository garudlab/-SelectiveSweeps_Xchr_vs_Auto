import msprime
import pyslim, tskit
import os.path
import numpy as np
import sys


def mutate_and_recap(inFile, outFile,r,N,mu,Q,hapsFile):
      # Load the .trees file from slim output
      ts = pyslim.load(inFile)
      #recapitate
      rts= pyslim.recapitate(ts,recombination_rate=r*Q, ancestral_Ne=int(N/Q))
      #simplify and sample
      individuals = rts.individuals_alive_at(0) #individuals alive at present
      sample = np.random.choice(individuals, size=100,replace=False)
      keep_nodes = []
      for i in sample:
           keep_nodes.extend(rts.individual(i).nodes)
      print("simplification")
      sts = rts.simplify(keep_nodes)

      print(f"Before, there were {rts.num_samples} sample nodes (and {rts.num_individuals} individuals)\n"
      f"in the tree sequence, and now there are {sts.num_samples} sample nodes\n"
      f"(and {sts.num_individuals} individuals).")

      #add mutations
      print("adding neutral mutations")
      mutated = pyslim.SlimTreeSequence(msprime.sim_mutations(sts,rate=mu*Q,model=msprime.SLiMMutationModel(type=1),keep=True,))


      print(f"The tree sequence now has {mutated.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {mutated.diversity()}, "
      f"and number of sites is {mutated.num_sites}.")

      ###output S and Pi ?
      
      
      print("output haplotypes")

      mutated2 = pyslim.convert_alleles(
            pyslim.generate_nucleotides(mutated))

      #get positions of variant sites
      pos=[]
      for variant in mutated2.variants():
            pos.append(int(variant.position))

      # get haplotypes
      haplo=[]
      for i, h in enumerate(mutated2.haplotypes(isolated_as_missing=None, missing_data_character="N", impute_missing_data=None)):
            haplo.append(list(h))

      haplo=np.array(haplo)
      n_rows = haplo.shape[0]
      haplo_100=haplo[np.random.choice(n_rows, size=100, replace=False), :]

      # get output with H12 format
      haplo_100= np.vstack([np.array(pos), haplo_100])
      out_haplo_100=np.transpose(haplo_100)
      np.savetxt(hapsFile, out_haplo_100, delimiter=",",fmt='%s')


      # print("output VCF")
      # try:
      #       #alive = mutated.individuals_alive_at(0)
      #       with open(outFile, "w") as vcffile:
      #             mutated2.write_vcf(vcffile)
      #             #mutated2.write_vcf(vcffile, individuals=alive)
      # except Exception as e:
      #       print ("Error:")
      #       print (e)


def main():

    inFile = sys.argv[1]
    outFile = sys.argv[2]
    r = float(sys.argv[3])
    N = int(sys.argv[4])
    mu = float(sys.argv[5])
    Q = float(sys.argv[6])
    hapsFile = sys.argv[7]

    mutate_and_recap(inFile, outFile,r,N,mu, Q,hapsFile) 

#run main
if __name__ == '__main__':
    main()