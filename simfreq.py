from cavities import freq_appoach,gw_binary_freq_signal,gw_binary_time_signal,heaviside,pi

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
f_max=fisco

if -log(M) < 5:
    f_max = 10**(log(M)+5)*fisco


l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
Rmax=rmax/l ; Rmin=rmin/l
B0 = 5.0 ; 
K=5; # number of radial modes

dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
hp,hx,t=gw_binary_time_signal(M,dist,f_min,deltaT)
tmax=t[-1] ; dt=t[1]-t[0]
deltaF = 1/tmax
TFhp,omega=gw_binary_freq_signal(M,dist,f_min,f_max,deltaF)

if -log(M) < 5:
    TFhp = TFhp*heaviside(-omega/2/pi+fisco,0.5)

#TFhp1=interp1d(omega,TFhp,kind='nearest')
#omega1=2*pi*linspace(0,f_max,10*len(TFhp))
print(TFhp,omega)
rmsP_f,TFdP_f,omega,dP_f,t_f=freq_appoach(TFhp,omega,cavity,Rmax,0)
print(rmsP_f)
plt.loglog(omega,abs(TFdP_f)**2)
filename='SIMFREQM{0}D{1}'.format(Mexp,den_f)+cavity+'.pkl'
with open(filename, 'wb') as f:
    dump([M,omega,TFdP_f,rmsP_f,t_f,dP_f], f,protocol=-1)
    f.close()

plt.show()