from cavities import RMS_compute, fchirp,cav_parameters
import numpy as np
from sys import argv

di=float(argv[1])
ds=float(argv[2])
num=int(argv[3])
Rmax=float(argv[4])
Rid=int(argv[5])
B0=5
l=1
cavity='TEM'
K=5
alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)
Mexp=np.loadtxt('mass.in')

Masse = 10**(-1*Mexp)

Deno = np.linspace(di,ds,num)
RMStemp=np.empty(len(Deno))
rep=np.empty(len(Mexp))
for j in range(len(Mexp)):
    fisco=2200/Masse[j]
    for i in range(len(Deno)):
        fmin=fisco/Deno[i]
        RMStemp[i]=RMS_compute(fchirp,fmin,fisco,cavity,Rmax,Masse[j])
    rep[j]=np.median(RMStemp)
title="results/simRay_{0:04d}.txt".format(Rid)
file = open(title, "w+")
np.savetxt(file,rep)
file.close()

