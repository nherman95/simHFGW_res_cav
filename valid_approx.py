from sympy.abc import t,x,y,z
from sympy import diff,init_printing
from sympy.utilities.lambdify import lambdify
import sympy
import numpy as np
from numpy import linspace
from math import log,factorial,pi
import matplotlib.pyplot as plt
from scipy.constants import c,G
from matplotlib import rcParams
rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 18, 8
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 20
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.5
rcParams['lines.markersize'] = 0.5
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 18

cavity='TM'
Mexp=5
den_f=25
M=10**(-Mexp); #Msun
print(M)
fisco=2200/(M); #hz
f_min = fisco/den_f;# start frequency of inspiral
MSUN_SI=1.9884099021470415e+30
PC_SI=3.085677581491367e+16
c=299792458
G=6.6743e-11
l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
Rmax=rmax/l ; Rmin=rmin/l
B0 = 5.0 ; 
K=5; # number of radial modes
dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
time=linspace(0,5.4269e-05,58271)
print(time[-1])
print(len(time))
mu=M/2*MSUN_SI*G/c**2
m=2*M*MSUN_SI*G/c**2
nu=1/4


init_printing()
taus=nu/5/m*t
phis=taus**(5/8)/nu
xs=1/4*taus**(-1/4)
hpan2s=-2*mu*xs*sympy.cos(2*phis)/dist/PC_SI
hxan2s=-2*mu*xs*sympy.sin(2*phis)/dist/PC_SI
hpan2s=diff(hpan2s,t,t)
hxan2s=diff(hxan2s,t,t)
print(hpan2s)
Pp=1/6*z**2*hpan2s
Px=1/6*z**2*hxan2s
Qp=1/3*z**2*hpan2s
Qx=1/3*z**2*hxan2s
nsum=10
for n in range (3,nsum+1):
    print(n)
    hpan2s=diff(hpan2s,t)
    Pp=Pp+(n-1)/factorial(n+1)*z**n*hpan2s
    Px=Px+(n-1)/factorial(n+1)*z**n*hxan2s
    Qp=Qp+(n)/factorial(n+1)*z**n*hpan2s
    Qx=Qx+(n)/factorial(n+1)*z**n*hxan2s
#hxx=Pp
hyy=-Pp
hxy=Px
hxz=-1/z*(x*Pp+y*Px)
hyz=-1/z*(x*Px-y*Pp)
hzz=1/z**2*((x**2-y**2)*Pp+2*x*y*Px)
htt=1/z**2*((x**2-y**2)*(2*Qp-Pp)+2*x*y*(2*Qx-Px))
#htx=-1/z*(x*Qp+y*Qx)
hty=-1/z*(x*Qx-y*Qp)
htz=1/z**2*((x**2-y**2)*(Qp)+2*x*y*(Qx))
trace=-htt+hzz
j1y=-diff(hyy,z)+diff(hyz,y)
j1z=-diff(hyz,z)+diff(hzz,y)
j2y=diff(htz,t)-diff(hxz,x)-diff(hyz,y)-diff(hzz,z)
j2z=-diff(hty,t)+diff(hxy,x)+diff(hyy,y)+diff(hyz,z)
j3y=diff(trace,z)/2
j3z=-diff(trace,y)/2
S1=lambdify([t,z],(diff(j1y,z)-diff(j1z,y)))
#print((diff(j1y,z)-diff(j1z,y)))
#print("one")
S=-diff(j1y+j2y+j3y,z)+diff(j1z+j2z+j3z,y)
Stota=lambdify([t,x,y,z],S)
print("SYMBOLIC DONE")
size=10
rarr=linspace(0,1,size)
tharr=linspace(0,2*pi,size)
rmesh,thmesh=np.meshgrid(rarr,tharr)
xmesh=rmesh*np.cos(thmesh)
ymesh=rmesh*np.sin(thmesh)
zp=[0.1,-0.5,1.75,-3]
k=50000
step=250


def Squad(zp):
    arr=np.array([])
    arr2=np.array([])
    for k in range(50000,200,-step):
        print(k)
        Smesh=np.zeros((size,size))
        for i in range(size):
            for j in range(size):
                Smesh[i,j]=Stota(c*time[k],xmesh[i,j],ymesh[i,j],zp)
        arr2=np.append(arr2,np.sqrt(pi)*abs(S1(c*time[k],zp)))
        arr=np.append(arr,np.sqrt(np.trapz(np.trapz(rmesh*Smesh[:,:]**2,tharr,axis=0),rarr,axis=0)))
        

    return arr,arr2


fig, ax=plt.subplots(1,len(zp))
fig2, ax2=plt.subplots(1,len(zp))
for i in range(len(zp)):
    print(i)
    Stot,S1quad=Squad(zp[i])
    ax[i].plot(Stot,label=r'$S$')
    ax[i].plot(7/6*S1quad,label=r"$\frac{7}{6}\,S_1$")
    ax2[i].plot(Stot/S1quad,label=r"$\frac{S}{S_1}$")
    ax2[i].plot(7/6*np.ones(len(Stot)),label=r"$\frac{7}{6}$")
    ax[i].set_ylim(0,1e-29)
    ax[i].set_xlabel(r'Time (arb. unit)')
    ax2[i].set_ylim(0,2)
    ax2[i].set_xlabel(r'Time (arb. unit)')
    #ax2[i].legend()
    ax[i].set_title(r"z={} m".format(zp[i]))
    ax2[i].set_title(r"z={} m".format(zp[i]))
ax[-1].legend(bbox_to_anchor =(1.45, 1.00))
ax[0].set_ylabel(r'$L^2$-norm of the source term ($T^{-1}m^{-2}$)')
ax2[-1].legend(bbox_to_anchor =(1.35, 1.00))
plt.subplots_adjust(wspace=0.25)
fig.savefig("approx.png")
fig2.savefig("approx2.png")
plt.show()