from cavities import RMS_compute, fchirp,gw_binary_freq_signal,gw_binary_time_signal,freq_appoach,time_approach
import numpy as np
from sys import argv
from math import log,pi

cavity='TM'
Mid=int(argv[1])
Mexp=float(argv[2])
di=float(argv[3])
ds=float(argv[4])
num=int(argv[5])

den_f=np.linspace(di,ds,num)
rmsP_f=np.empty(num)
rmsP_t=np.empty(num)
rmsP_rir=np.empty(num)

for i in range(len(den_f)):
    M=10**(-Mexp); #Msun
    fisco=2200/(M); #hz
    f_min = fisco/den_f[i];# start frequency of inspiral
    f_max=fisco

    if -log(M) < 5:
        f_max = 10**(log(M)+5)*fisco


    l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
    Rmax=rmax/l ; Rmin=rmin/l
    B0 = 5.0 ; 
    K=5; # number of radial modes

    dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
    hp,_,t=gw_binary_time_signal(M,dist,f_min,deltaT)
    tmax=t[-1] ; dt=t[1]-t[0]
    deltaF = 1/tmax
    TFhp,omega=gw_binary_freq_signal(M,dist,f_min,f_max,deltaF)

    if -log(M) < 5:
        TFhp = TFhp*np.heaviside(-omega/2/pi+fisco,0.5)

    rmsP_t[i],_=time_approach(hp,t,cavity,Rmax)
    rmsP_f[i],_,_,_=freq_appoach(TFhp,omega,cavity,Rmax)
    rmsP_rir[i]=RMS_compute(fchirp,f_min,fisco,cavity,Rmax,M)
rep=np.array(np.array([np.median(rmsP_t),np.median(rmsP_f),np.median(rmsP_rir)]))
title="results/simMass_{0:03d}.txt".format(Mid)
file = open(title, "w+")
np.savetxt(file,np.matrix(rep))

