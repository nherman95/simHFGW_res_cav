from cavities import Sh_toy, Sh_toyff1, Sh_toyfm1,RMS_compute,freq_appoach,SSh_toy,SSh_toyff1,SSh_toyfm1
import matplotlib.pyplot as plt
from numpy import pi,linspace,sqrt
from scipy.stats import linregress
import numpy as np
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter,NullFormatter
rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 9.5, 8
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 36
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.5
rcParams['lines.markersize'] = 0.5
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22

l=1
cavity='TEM'
Rmax=5
K=5
omegaGW=1e-10
fmax =1e10
fmin = 1e5
dt = 0.25*fmax**-1
nbrpt = 1/fmin/dt
freq = linspace(fmin,fmax,int(nbrpt))
omega = 2*pi*freq
fcutarr=[1e7,1e8,1e9]
psdarr=[]

for fcut in fcutarr:
     test=sqrt(2)*SSh_toy(freq,omegaGW,fcut)
     rmsP_f,TFdP_f_p,omega,dP_f,tdE_f=freq_appoach(test,omega,cavity,Rmax)
     psdarr.append(abs(TFdP_f_p)**2)
#plt.loglog(freq,abs(TFdP_f_p)**2,'b')
plt.loglog(freq,psdarr[0],'--r',label=r'$\nu_{cut}=10^7$ Hz')
plt.loglog(freq,psdarr[1],'r',linewidth=5,label=r'$\nu_{cut}=10^8$ Hz')
plt.loglog(freq,psdarr[2],':r',label=r'$\nu_{cut}=10^9$ Hz')
plt.xlim(1e5,1e10)
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD of induced Power ($W^2$/Hz)')
plt.legend(loc=2)
plt.savefig("psd.png")
# print(RMS_compute(Sh_toy,fmin,fmax,cavity,Rmax))
# plt.figure()
# # plt.loglog(freq,Infl,'r')#,label='inflation')
# # #plt.loglog(freq,abs(TFdP_f_p),'r')#,label='inflation')
# # #plt.loglog(freq,4*pi*fct_bode(freq,cavity,Rmax)*Sh_infl(freq))
# # plt.loglog(freq,Phasetran,'g',label='phase transition')
# # plt.loglog(freq,preheat,'b',label='pre-heating')
# plt.loglog(freq,Sh_toy(freq),'b',label='pre-heating')
# #plt.loglog(freq,fchirp(freq,1e-5))
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('$h_c$')
# #plt.legend()
# plt.savefig("psd.png")

# plt.figure()
# rmsP_f,TFdP_f,dP_f,t_f=freq_appoach(Infl,omega,fmin, fmax,cavity,Rmax)
# plt.loglog(freq,abs(TFdP_f),'r')
# rmsP_f,TFdP_f,dP_f,t_f=freq_appoach(Phasetran,omega,fmin, fmax,cavity,Rmax)
# plt.loglog(freq,abs(TFdP_f),'g')
# rmsP_f,TFdP_f,dP_f,t_f=freq_appoach(preheat,omega,fmin, fmax,cavity,Rmax)
# plt.loglog(freq,abs(TFdP_f),'b')
# plt.show()


#rmaxarr = 10**(np.linspace(0.101,1,5))
#rmaxarr = linspace(1,10,100)
rmaxarr = 10**(np.array(list(linspace (0,0.2,5)) + list(linspace (0.2,1,7))))
print(rmaxarr)
fmin=1
fmax=1e9
RMS_toy=[]
RMS_toy2=[]
RMS_toy3=[]
#RMS_toy4=[]
#RMS_toy5=[]
fcut=1e8
# RMS_pt=[]
# RMS_ph=[]
for r in range(len(rmaxarr)):
     print(r)
     l = 1.0  ;rmax = rmaxarr[r] ; rmin = 0.1 #m
     Rmax=rmax/l ; Rmin=rmin/l
     RMS_toy.append(RMS_compute(SSh_toy,fmin,fmax,cavity,Rmax,omegaGW,fcut))
     RMS_toy2.append(RMS_compute(SSh_toyfm1,fmin,fmax,cavity,Rmax,omegaGW,fcut))
     RMS_toy3.append(RMS_compute(SSh_toyff1,fmin,fmax,cavity,Rmax,omegaGW,fcut))
     #RMS_toy4.append(RMS_compute(Sh_toyfm2,fmin,fmax,cavity,Rmax,omegaGW,fcut))
     #RMS_toy5.append(RMS_compute(Sh_toyff2,fmin,fmax,cavity,Rmax,omegaGW,fcut))
#     RMS_pt.append(RMS_compute(Sh_pt,fmin,fmax,cavity,Rmax))
#     RMS_ph.append(RMS_compute(Sh_ph,fmin,fmax,cavity,Rmax))


    
# #print(time()-t0)
# ## HFSGWB ps des modèles
# ##Plot SGWB en fct de rmax
# #regime linéaire avant fcut
# slope_infl, intercept_infl, r_value_infl, p_value_infl, std_err_infl = linregress(log10(rmaxarr[70:]), log10(RMS_toy[70:]))
# # slope_pt, intercept_pt, r_value_pt, p_value_pt, std_err_pt = linregress(log10(rmaxarr[70:]), log10(RMS_pt[70:]))
# # slope_ph, intercept_ph, r_value_ph, p_value_ph, std_err_ph = linregress(log10(rmaxarr[70:]), log10(RMS_ph[70:]))
# # print(r_value_infl)
# # print(r_value_pt)
# # print(r_value_ph)
# # plt.loglog(rmaxarr,rmaxarr**(slope_infl)*10**(intercept_infl),linestyle='--',color='r')
# # plt.loglog(rmaxarr,rmaxarr**(slope_pt)*10**(intercept_pt),linestyle='--',color='g')
# # plt.loglog(rmaxarr,rmaxarr**(slope_ph)*10**(intercept_ph),linestyle='--',color='b')
# #print(RMS_pt)
plt.figure()
plt.loglog(rmaxarr,RMS_toy,'^',markersize=10,color='r',label=r'$\Omega_{GW}$ constant')
plt.loglog(rmaxarr,RMS_toy2,'^',markersize=10,color='g',label=r'$\Omega_{GW} \propto \nu^{-1}$ ')
plt.loglog(rmaxarr,RMS_toy3,'^',markersize=10,color='b',label=r'$\Omega_{GW} \propto \nu$')
#plt.loglog(rmaxarr,RMS_toy4,'^',markersize=10,color='y')
#plt.loglog(rmaxarr,RMS_toy5,'^',markersize=10,color='k')
# plt.loglog(rmaxarr,RMS_pt,'.',markersize=3,color='g',label='phase transition')
# plt.loglog(rmaxarr,RMS_ph,'.',markersize=3,color='b',label='pre-heating')
slope_infl, intercept_infl, r_value_infl, p_value_infl, std_err_infl = linregress(np.log10(rmaxarr[-5:]), np.log10(RMS_toy[-5:]))
print(slope_infl,r_value_infl)
slope_infl, intercept_infl, r_value_infl, p_value_infl, std_err_infl = linregress(np.log10(rmaxarr[-5:]), np.log10(RMS_toy2[-5:]))
print(slope_infl,r_value_infl)
slope_infl, intercept_infl, r_value_infl, p_value_infl, std_err_infl = linregress(np.log10(rmaxarr[-5:]), np.log10(RMS_toy3[-5:]))
print(slope_infl,r_value_infl)
#slope_infl, intercept_infl, r_value_infl, p_value_infl, std_err_infl = linregress(np.log10(rmaxarr[-5:]), np.log10(RMS_toy4[-5:]))
#print(slope_infl,r_value_infl)
#slope_infl, intercept_infl, r_value_infl, p_value_infl, std_err_infl = linregress(np.log10(rmaxarr[-5:]), np.log10(RMS_toy5[-5:]))
#print(slope_infl,r_value_infl)
plt.xlabel('Outer radius (m)')
plt.ylabel('Induced RMS power (W)')
plt.xlim(1,10)
ax=plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(NullFormatter())
#plt.xticks([])
plt.xticks([1,5,10],['1','5','10'])
#plt.yticks()
plt.legend()
plt.savefig("radius.png")
