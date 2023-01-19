from numpy import loadtxt,arange,pi,sqrt
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 10,10
rcParams['axes.labelsize'] = 28
rcParams['axes.titlesize'] = 36
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 8
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
Mexp=loadtxt('mass.in')
M=10**(-1*Mexp[6:])
Y=loadtxt('simMass.txt')
fig=plt.figure()
ax=plt.gca()
#Y=Y*4.640850738526967e-26*72511535.92617007*2*pi*5
plt.loglog(M,Y[:,0],'b',label="time domain")
plt.loglog(M,Y[:,1],'r',label="frequency domain")
plt.loglog(M,Y[:,2],'g',label="RMS power integral estimation")
plt.loglog(M[:int(len(M)/3)],1e-21*M[:int(len(M)/3)]**(-11/6),'--k',linewidth=1)
plt.annotate(r"$\propto M^{-11/6}$",(3e-4,3e-15), xytext=(3e-5,1e-16),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=15,fontweight='bold')
plt.annotate(r"$\propto M^{7/6}$",(1e-7,7e-16), xytext=(2e-8,1e-14),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=15,fontweight='bold')
plt.loglog(M[2*int(len(M)/3):],1e-7*M[2*int(len(M)/3):]**(7/6),':k',linewidth=1)
# plt.fill_between(M,Y[1,:],Y[2,:],'b',alpha=0.35)
plt.axis([1e-8,1e-3,1e-17,1e-9])
xt=10.**arange(-8,-2)
print(xt)
plt.xticks(xt)
ax.set_box_aspect(1)
plt.xlabel(r"Mass ($M_\odot$)")
plt.ylabel(r"RMS Induced Power (W)")
plt.title(r'TM Cavity',fontweight='bold')
plt.legend()
plt.savefig("pbhmass.pdf")

nii=loadtxt('niikura.csv')
plt.figure()
ax=plt.gca()
test=(1e-14*4.2*M**(0.29)/Y[:,2]/1000)**(3/2)
test2=(1e-14*58/Y[:,2]/1000)**(3/2)
plt.loglog(M,test,label='Primordial binaries')
plt.fill_between(M,test,y2=1e1,alpha=0.3)
plt.loglog(M,test2,label='Tidal capture in clusters')
plt.fill_between(M,test2,y2=1e1,alpha=0.5)
plt.loglog(nii[:,0],nii[:,1],marker='.',alpha=0.7,label="OGLE+Gaia: Niikura 19'")
plt.ylim(1e-12,1)
plt.xlim(1e-8,1e-3)
ax.set_box_aspect(1)
plt.ylabel(r"$\tilde{f}_{\rm PBH}$")
plt.xlabel(r"Mass ($M_\odot$)")
plt.legend()
plt.savefig("niikura.pdf")