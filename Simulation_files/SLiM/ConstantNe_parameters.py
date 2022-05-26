import numpy as np
import sys
#from scipy.stats import rv_discrete
import os.path

rho_in=float(sys.argv[1])
adaptive_theta=float(sys.argv[2])
selection=sys.argv[3] # boolean indicating if we are similating with selection
Ne=float(sys.argv[4]) 
Q=int(sys.argv[5])
chrType= sys.argv[6]
program=sys.argv[7]
run=sys.argv[8]
ite=sys.argv[9]

sb=np.random.uniform(0,1)
PF=np.random.uniform(0,1)
H=np.random.uniform(0,1)
AgeSweep=int(np.random.uniform(0,0.001)*4*Ne/Q) #multiply*0.75 for X



if selection == 'True':
    command='slim ' + ' -d R='+ str(rho_in) + ' -d N=' + str(Ne) + ' -d THETA_A=' +str(adaptive_theta) +' -d Q='+ str(Q) +' -d "ChrType=' +"'"+ str(chrType)+"'"+str('"')
    command += ' -d sb='+str(sb)+ ' -d PF='+str(PF)+ ' -d H='+str(H)+' -d AgeSweep='+str(AgeSweep)+ ' -d run='+str(run)+str(ite)+' ConstantNe_rescaling.slim'
    print(command)
else:
    sb=0
    command='slim ' + ' -d R='+ str(rho_in) + ' -d N=' + str(Ne) + ' -d THETA_A=' +str(adaptive_theta) +' -d Q='+ str(Q) +' -d "ChrType=' +"'"+ str(chrType)+"'"+str('"')
    command += ' -d sb='+str(sb)+ ' -d PF='+str(PF)+ ' -d AgeSweep='+str(AgeSweep)+ ' -d run='+str(run)+str(ite)+' ConstantNe_rescaling.slim'
    print(command)


