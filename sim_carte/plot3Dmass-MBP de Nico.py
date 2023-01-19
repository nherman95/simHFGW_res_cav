from numpy import loadtxt,array,log10,arange,meshgrid,empty
import matplotlib.pyplot as plt
from matplotlib import rcParams,colors
from numpy.core.function_base import linspace
rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 12, 10
rcParams['axes.labelsize'] = 28
rcParams['axes.titlesize'] = 36
rcParams['font.size'] = 12
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 8
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14
Mexp=loadtxt('mass.in')
R=loadtxt('ray.in')
Y=loadtxt('simtotal.txt')
fig=plt.figure()
ax=plt.gca()
x = []
y = []
z = []
for m in range (len(Mexp)):
    for r in range (len(R)):
        x.append(10**(-Mexp[m]))
        y.append((R[r]))
        z.append(Y[m+r*len(R)]) #abs(save_A_fit[m,r]/Y[m,r])
        
x = array(x)
y = array(y)
z = array(z)
plt.scatter(x,y,marker='s',c = z, cmap='jet',norm=colors.LogNorm(1e-20,1e-10))#seismic#coolwarm
cbar=plt.colorbar()
cbar.ax.set_title("$P_{RMS}$ (W)",fontsize=16)
cbar.ax.tick_params(labelsize=12)

cbar.set_ticks(10.0**linspace(-20,-10,11))
plt.xscale('log')
plt.xlim(10**-8,10**-2)
plt.xticks(10.0**linspace(-8,-2,7))
plt.yticks(linspace(1,10,10))
plt.ylim(1,10)
plt.xlabel('Mass($M_\odot$)')
plt.ylabel('Cavity radius (m)')
#M,r=meshgrid(10**(-Mexp),R)
#plt.tricontour(x,y,z,[1e-14,1e-11,1e-10,1e-9],colors='k')
ax.set_box_aspect(1)
plt.savefig("carte.pdf")
