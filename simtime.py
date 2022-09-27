from cavities import gw_binary_time_signal,time_approach
from numpy import linspace
from math import log
from pickle import dump
import matplotlib.pyplot as plt
cavity='TM'
Mexp=3
den_f=25
M=10**(-Mexp); #Msun
print(M)
fisco=2200/(M); #hz
f_min = fisco/den_f;# start frequency of inspiral


l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
Rmax=rmax/l ; Rmin=rmin/l
B0 = 5.0 ; 
K=5; # number of radial modes

dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
hp,hx,t=gw_binary_time_signal(M,dist,f_min,deltaT)
tmax=t[-1] ; dt=t[1]-t[0]

rmsP_t,dP_t=time_approach(hp,t,cavity,Rmax)
t=linspace(0,t[-1],len(dP_t))
plt.plot(t,dP_t)
print(rmsP_t)
plt.show()

filename='SIMTIMEM{0}D{1}'.format(Mexp,den_f)+cavity+'.pkl'
with open(filename, 'wb') as f:
    dump([M,t,dP_t,rmsP_t], f,protocol=-1)
    f.close()