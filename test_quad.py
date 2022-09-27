from math import log2,pi
from cavities import *
import numpy as np
from scipy.integrate import quadrature,romberg,simps,romb,trapz
from scipy.stats import norm
from scipy.special import jv
import matplotlib.pyplot as plt
# M=1e-5
# npts=1e7
# fmax=2200/M
# fmin=1
# print(fmax-fmin)
# #print(sqrt(quadrature(f_RMS_mass_deno,fmin,fmax,args=(M,'TM',5),tol=1e-15,rtol=1e-15)))
# #print(sqrt(romberg(f_RMS_mass_deno,fmin,fmax,args=(M,'TM',5),tol=1e-15,rtol=1e-15)))
# #fspan=linspace(fmin,fmax,int(npts))
# fspan2=linspace(fmin,fmax,2**round(log2(npts))+1)
# #print(sqrt(simps(f_RMS_mass_deno(fspan,M,'TM',5),fspan)))
# #print(sqrt(trapz(f_RMS_mass_deno(fspan,M,'TM',5),fspan)))
# # test=norm.pdf(fspan2,fspan2[2**10],1e-12)
# # plt.plot(test)
# # plt.show()
# # print(sqrt(romb(test**2,fspan2[1]-fspan2[0]))/(fmax-fmin))
# #print(round(log2(npts)))
# plt.loglog(fspan2,fchirp(fspan2,M))
# plt.show()
# print(RMS_compute(fchirp,fmin,fmax,'TM',5,M))
# #plt.loglog(fspan,f_RMS_mass_deno(fspan,M,'TM',5))
# #plt.loglog(fspan2,f_RMS_mass_deno(fspan2,M,'TM',5))
# #plt.show()

Rmax=5
alphak,Ir,Ith,normr,normth=cav_parameters('TM',5,1,Rmax)
rlin=np.linspace(1e-15,Rmax,1000)
for K in range(5):
    f_nr=(jv(1,alphak[K]*rlin)/rlin)
    nr=trapz(f_nr,rlin)
    print(nr)
    print(abs(Ir[K]-nr))
