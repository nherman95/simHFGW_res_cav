from cavities import fct_bode,cav_parameters,time_approach,RMS_compute,fchirp
import matplotlib.pyplot as plt
import numpy as np

pi=np.pi
B0=5
l=1
cavity='TEM'
Rmax=5
K=5
#alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)

fplot=np.logspace(5,9,1000000)
#print(fplot)
plt.loglog(fplot,fct_bode(fplot,cavity,Rmax))
plt.savefig("bodeTEM.png")
# test=[]
# spot
# # for i in range(len(fplot)):
# #     dt=1/20/fplot[i]
# #     print(fplot[i])
# #     t = np.linspace(0,10/fplot[i],round(10/fplot[i]/dt))
        
# #     tmax = t[-1]
# #     hp = np.cos(2*pi*fplot[i]*t)+np.sin(2*pi*fplot[i]*t)#*4.640850738526967e-31*(72511535.92617007/f)
# #     hGW=1; hp=hp/hGW;
# #     rmsP_t,dP_t=time_approach(hp,t,cavity,Rmax)
# #     test.append(rmsP_t)

# # plt.loglog(fplot,fct_bode(fplot,cavity,Rmax))
# # plt.figure()
# # plt.plot(fplot,fct_bode(fplot,cavity,Rmax))
# # plt.loglog(fplot,test)
# # plt.show()

# Mexp=np.linspace(1,9,30)

# Masse = 10**(-1*Mexp)

# #Deno = np.linspace(di,ds,num)
# #RMStemp=np.empty(len(Deno))
# rep=np.empty(len(Mexp))
# for j in range(len(Mexp)):
#     print(j)
#     fisco=2200/Masse[j]
#     fmin=fisco/25
#     rep[j]=RMS_compute(fchirp,fmin,fisco,cavity,Rmax,Masse[j])

# plt.loglog(Masse,rep)
# plt.show()