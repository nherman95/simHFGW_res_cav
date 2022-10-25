from cavities import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.ticker import ScalarFormatter,NullFormatter
from matplotlib.patches import Ellipse
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
rcParams['lines.linewidth'] = 3.0
rcParams['lines.markersize'] = 8
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22

pi=np.pi
B0=5
l=1
cavity='TM'
Rmax=5
K=5
alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)

fplot=np.linspace(1e5,1e9,2000000)
plt.loglog(fplot,fct_bode(fplot,cavity,Rmax),'b')
plt.xlim(1e5,1e9)
plt.xlabel('Frequency (Hz)')
plt.ylabel('$P_{RIR}$ (W)')
plt.savefig("bodeTM.png")
cavity='TEM'
plt.figure()
print(fct_bode(fplot,cavity,Rmax))
plt.loglog(fplot,fct_bode(fplot,cavity,Rmax),'b')
plt.xlim(1e5,1e9)
plt.xlabel('Frequency (Hz)')
plt.ylabel('$P_{RIR}$ (W)')
plt.savefig("bodeTEM.png")

plt.figure()
prmslim=1e-14
plt.loglog(fplot,sqrt(fplot)*Sh_toy(fplot,1e-10,1e8),'b',label='stochastic sources')
chirp=[]
for i in range(len(fplot)):
        chirp.append(fchirp(fplot[i],2200/fplot[i]))
        #chirp.append(fchirp(fplot[i],1e-5))

plt.loglog(fplot,2*fplot*array(chirp),'r',label='inspiral sources')
#plt.loglog(fplot,sqrt(fplot)*prmslim/fct_bode(fplot,cavity,Rmax),label='stochastic')
plt.loglog(fplot,prmslim/fct_bode(fplot,cavity,Rmax),':k')
plt.legend()
plt.xlabel('Frequency (Hz)')
plt.ylabel('Characteristic Strain')
rng = plt.axis()
x_scale = 9.5 * 0.78 / (rng[1] - rng[0])         
y_scale = 8 * 0.80 / (rng[3] - rng[2])  
# angle= np.degrees(np.arctan2((np.log(prmslim/fct_bode(fplot[1],cavity,Rmax))-np.log(prmslim/fct_bode(fplot[0],cavity,Rmax))),np.log(fplot[1])-np.log(fplot[0])))
# print(angle)
l=np.array([4e5,1e-30])
# transangle=plt.gca().transData.transform_angles(np.array((angle,)),l.reshape((1, 2)))[0]
# print(transangle)
plt.annotate("Detectable strain at $\mathbf{10^{-14}}$ W ",l,rotation=-61.5/1.75, fontweight='bold',fontsize=18)
plt.xlim(1e5,1e9)
ell=Ellipse((9e7,1e-33),1.5e8,1e-30,alpha=0.35,color='gold')
ax=plt.gca()
ax.add_artist(ell)
l=np.array([1.5e6,1e-33])
plt.annotate("Resonant modes",l,color='goldenrod',fontweight='bold',fontsize=18)
plt.savefig("detlim.png")


# Masse = np.linspace(2,8,50); Masse = 10**(-1*Masse)

# Deno = np.array([25])#np.linspace(10,35,5)
# RMSlist=np.empty(len(Masse))
# RMSlistQ1=np.empty(len(Masse))
# RMSlistQ3=np.empty(len(Masse))
# for j in range(len(Masse)):
#         print(j)
#         fisco=2200/Masse[j]
#         RMStemp=np.empty(len(Deno))
#         for i in range(len(Deno)):
#                 fmin=fisco/Deno[i]
#                 RMStemp=RMS_compute(fchirp,fmin,fisco,'TM',5,Masse[j])
#         RMSlist[j]=np.median(RMStemp)
#         RMSlistQ1[j]=np.percentile(RMStemp,25)
#         RMSlistQ3[j]=np.percentile(RMStemp,75)


        
        
#print(RMStest)
#for k in range (np.shape(RMStest)[1]):
#    RMStest[:,k]=np.sort(RMStest[:,k])
#bi=round(10/100*len(RMStest[:,0]))
#bs=round(90/100*len(RMStest[:,0]))
# RMSlist=[]
# for k in range (len(Masse)):
#     if 10**-4.9<Masse[k]<10**-4:
#         RMSlist.append(mean(RMStest[0:int(600/1000*len(RMStest[:,0])*k**2/1400),k]))
#     else :
#         RMSlist.append(mean(RMStest[0:int(900/1000*len(RMStest[:,0])),k]))
#RMSlist=np.mean(RMStest,0)
#RMSlist2=np.median(RMStest,0)
#RMSlist2Q1=np.percentile(RMStest,25,0)
#RMSlist2Q3=np.percentile(RMStest,75,0)
#print(RMSlist)
# plt.figure()
# plt.loglog(Masse,RMSlist)
# #plt.loglog(Masse,RMSlist2)
# plt.loglog(Masse,RMSlistQ1)
# plt.loglog(Masse,RMSlistQ3)
# plt.savefig("masse.png")
# #plt.show()
