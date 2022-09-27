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
rcParams['figure.figsize'] = 15.5, 7
rcParams['axes.labelsize'] = 28
rcParams['axes.titlesize'] = 36
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 8
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
Mexp=loadtxt('mass.in')
M=10**(-1*Mexp)
Y=loadtxt('results/sim.txt')
fig=plt.figure()
ax=plt.gca()
#Y=Y*4.640850738526967e-26*72511535.92617007*2*pi*5
plt.loglog(M,Y[0,:],'b',linewidth=3)

plt.loglog(M,Y[1,:],'b--')
plt.loglog(M,Y[2,:],'b--')
plt.fill_between(M,Y[1,:],Y[2,:],'b',alpha=0.35)
plt.axis([1e-8,1e-2,1e-17,1e-9])
xt=10.**arange(-8,-1)
print(xt)
plt.xticks(xt)
ax.set_box_aspect(1)
plt.xlabel("Mass ($M_\odot$)")
plt.ylabel("RMS Induced Power (W)")
plt.title('TM Cavity',fontdict={'fontweight':'bold'})
plt.savefig("masse.png")
#plt.show()