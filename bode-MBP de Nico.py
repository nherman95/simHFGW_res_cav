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
#rcParams['text.usetex'] = True

pi=np.pi
B0=5
l=1
cavity='TM'
Rmax=5
K=5
alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)

fplot=np.linspace(1e5,1e9,2000000)
plt.loglog(fplot,fct_bode(fplot,cavity,Rmax),'b',label='TM')
plt.annotate("010 mode",(c*alphak[0]/2/pi,fct_bode(c*alphak[0]/2/pi,cavity,Rmax)*1.05),xytext=(5e6,1e22),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=15)
for k in range(1,len(alphak)):
        plt.annotate("{0}10".format(k),(c*alphak[k]/2/pi*0.9,fct_bode(c*alphak[k]/2/pi,cavity,Rmax)*1.2),fontsize=8)
plt.xlim(1e5,1e9)
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('$P_{RIR}$ (W)')
# plt.savefig("bodeTM.png")
cavity='TEM'
#plt.figure()
#print(fct_bode(fplot,cavity,Rmax))
plt.loglog(fplot,fct_bode(fplot,cavity,Rmax),'r',label='TEM')
plt.loglog(fplot[:int(len(fplot)/100)],2e-5*fplot[:int(len(fplot)/100)]**3,'--k',linewidth=1)
plt.annotate(r" $\propto \omega^3$ ",(1e6,2e13), xytext=(2e5,1e14),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=15)
plt.annotate(r" $\sin\left(\frac{\omega c}{2L}\right)$ ",(2e8,1e14), xytext=(3e7,1e13),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=15)

plt.xlim(1e5,1e9)
plt.xlabel('Frequency (Hz)')
plt.ylabel('$P_{RIR}$ ($W\cdot$ strain)')
plt.legend()
#plt.show()
plt.savefig("bodeTETEM.pdf")

plt.figure()
prmslim=1e-14
plt.loglog(fplot,sqrt(fplot)*SSh_toy(fplot,1e-10,1e8),'b',label='stochastic sources')
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
l=np.array([6e5,1e-30])
# transangle=plt.gca().transData.transform_angles(np.array((angle,)),l.reshape((1, 2)))[0]
# print(transangle)
plt.annotate("Detectable strain at $\mathbf{10^{-14}}$ W ",l,rotation=-61.5/1.75, fontweight='bold',fontsize=18)
plt.xlim(1e5,1e9)
ell=Ellipse((9e7,1e-33),1.5e8,1e-30,alpha=0.35,color='gold')
ax=plt.gca()
ax.add_artist(ell)
l=np.array([1.5e6,1e-33])
plt.ylim(1e-38,1e-24)
plt.annotate("Resonant modes",l,color='goldenrod',fontweight='bold',fontsize=18)
plt.savefig("detlim.pdf")
plt.figure()
plt.loglog(fplot,sqrt(fplot)*SSh_toy(fplot,1e-10,1e8),'b')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Characteristic Strain')
plt.xlim(1e5,1e9)
plt.ylim(1e-33,1e-28)
plt.savefig("toymodel.pdf")

fig, ax=plt.subplots(1,K,figsize=(18,7))
fig2, ax2=plt.subplots(1,K,figsize=(18,7))
r=linspace(1e-16,5,200)
phi=linspace(-pi,pi,200)
R,P=np.meshgrid(r,phi)
X=R*cos(P)
Y=R*sin(P)
Z=np.empty(np.shape(R))
Z2=np.empty(np.shape(R))
print(np.shape(R))
for k in range(K):
        print(k)
        for i in range(np.size(R[0,:])):
                Z[i,:]=jv(1,alphak[k]*R[i,:])/R[i,:]*cos(P[i,:])/normk[k]
                Z2[i,:]=-(alphak[k]*jv(0,alphak[k]*R[i,:])-jv(1,alphak[k]*R[i,:])/R[i,:])*sin(P[i,:])/normk[k]
        cp=ax[k].scatter(X,Y,c=Z/np.max(np.abs(Z)),cmap='jet')
        cp.set_clim(-1,1)
        ax[k].set_box_aspect(1)
        ax[k].set_xticks([])
        ax[k].set_yticks([])
        ax[k].set_title("k={0:d}".format(k))
        cp2=ax2[k].scatter(X,Y,c=Z2/np.max(np.abs(Z2)),cmap='jet')
        cp2.set_clim(-1,1)
        ax2[k].set_box_aspect(1)
        ax2[k].set_xticks([])
        ax2[k].set_yticks([])
        ax2[k].set_title("k={0:d}".format(k))
        #cp2=ax2[k].scatter(X,Y,c=Z2/np.max(np.abs(Z2)),cmap='jet')
        #cp2.set_clim(-1,1)
        #ax2[k].set_box_aspect(1)
fig.colorbar(cp,ax=ax,location='bottom')  
fig.suptitle("Representation of $\psi^r_{k10}$",fontsize=48)  
fig2.colorbar(cp2,ax=ax2,location='bottom')  
fig2.suptitle("Representation of $\psi^\\theta_{k10}$",fontsize=48)  
fig.savefig("psir.pdf")
fig2.savefig("psith.pdf")



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
