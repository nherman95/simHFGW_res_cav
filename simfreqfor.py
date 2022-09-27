from numpy import append,mean,var
from cavities2 import *
import matplotlib.pyplot as plt
cavity='TM'
Mexp=linspace(2,8,51)
den_f=linspace(24,26,10)
RMS_mean=array([])
RMS_ste=array([])
for i in range(len(Mexp)):
    M=10**(-Mexp[i]); #Msun
    print(M)
    fisco=2200/(M); #hz
    Stock=[]
    for j in range(len(den_f)):
        f_min = fisco/den_f[j];# start frequency of inspiral
        f_max=fisco

        if -log(M) < 5:
            f_max = 10**(log(M)+5)*fisco


        l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
        Rmax=rmax/l ; Rmin=rmin/l
        B0 = 5.0 ; 
        K=5; # number of radial modes

        dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
        hp,t=gw_binary_time_signal(M,dist,f_min,deltaT)
        tmax=t[-1] ; dt=t[1]-t[0]
        deltaF = 1/tmax
        TFhp,omega=gw_binary_freq_signal(M,dist,f_min,f_max,deltaF)

        if -log(M) < 5:
            TFhp = TFhp*heaviside(-omega/2/pi+fisco,0.5)

        rmsP_f,TFdP_f,dP_f,t_f=freq_appoach(TFhp,omega,f_min,f_max,cavity,Rmax)
        Stock.append(rmsP_f)
    RMS_mean=append(RMS_mean,mean(Stock))
    RMS_ste=append(RMS_ste,sqrt(var(Stock)/len(Stock)))

plt.loglog(10**-Mexp,RMS_mean,color='b')
# plt.loglog(10**-Mexp,RMS2_mean,color='r')
print(len(RMS_mean))
print(len(RMS_ste))
#plt.loglog(10**-Mexp,RMS_mean+RMS_ste,color='b',linestyle='--')
#plt.loglog(10**-Mexp,RMS_mean-RMS_ste,color='b',linestyle='--')
# plt.loglog(10**-Mexp,RMS2_mean+RMS2_var,color='r',linestyle='--')
# plt.loglog(10**-Mexp,RMS2_mean-RMS2_var,color='r',linestyle='--')

plt.loglog(10**-Mexp[int(len(Mexp)/1.6):],10**(2*Mexp[int(len(Mexp)/1.6):])*3.3/10**20,'--',label='$M^{-2}$',color='grey')
plt.loglog(10**-Mexp[0:int(len(Mexp)/2.5)],10**-Mexp[0:int(len(Mexp)/2.5)]*0.28/10**5,'--',label='$M^{-1}$',color='grey')
plt.legend()
plt.savefig("test.png")


        