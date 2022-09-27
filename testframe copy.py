from cavities import *
from sympy.abc import t,x,y,z
from sympy import diff,init_printing,simplify
from sympy.utilities.lambdify import lambdify
import sympy
import numpy as np
from numpy import linspace,flip
from math import log,factorial
from pickle import dump
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
from matplotlib.animation import *
cavity='TM'
Mexp=5
den_f=25
M=10**(-Mexp); #Msun
print(M)
fisco=2200/(M); #hz
f_min = fisco/den_f;# start frequency of inspiral


l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
Rmax=rmax/l ; Rmin=rmin/l
B0 = 5.0 ; 
K=5; # number of radial modes

dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
hp,hx,time=gw_binary_time_signal(M,dist,f_min,deltaT)
Mc=M*MSUN_SI**(6/5)/(2*M*MSUN_SI)**(1/5)
#
tau0=2.18*(1.21*MSUN_SI/Mc)**(5/8)*(100/f_min)**(8/3)
#tauf=2.18*(1.21*MSUN_SI/Mc)**(5/8)*(100/fisco)**(8/3)
#t#au=flip(linspace(t[1],tau0,len(t)))
# tau=flip(t[1:])
# print(tau)
# phi=-2*(5*G*Mc/c**3)**(-5/8)*tau**(5/8)
# hpan=1/dist/PC_SI*(G*Mc/c**2)**(5/4)*(5/c/tau)**(1/4)*cos(phi)

mu=M/2*MSUN_SI*G/c**2
m=2*M*MSUN_SI*G/c**2
nu=1/4
tau=c*nu/5/m*flip(time[1:])
phi=tau**(5/8)/nu
xx=1/4*tau**(-1/4)
hpan2=-2*mu*xx*np.cos(2*phi)/dist/PC_SI
#plt.plot(time[:-1]*c,hpan2)
#plt.show()

init_printing()
taus=nu/5/m*t
phis=taus**(5/8)/nu
xs=1/4*taus**(-1/4)
hpan2s=-2*mu*xs*sympy.cos(2*phis)/dist/PC_SI
hxan2s=-2*mu*xs*sympy.sin(2*phis)/dist/PC_SI
#hpan2s=sympy.cos(t**(5/8))
#hxan2s=sympy.sin(t**(5/8))
hpan2s=diff(hpan2s,t,t)
hxan2s=diff(hxan2s,t,t)
print(hpan2s)
Pp=1/6*z**2*hpan2s
Px=1/6*z**2*hxan2s
Qp=1/3*z**2*hpan2s
Qx=1/3*z**2*hxan2s
nsum=5
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
#dondert=lambdify([t,x,y,z],-diff(htt,t)+diff(htx,x)+diff(hty,y)+diff(htz,z))
#donderx=lambdify([t,x,y,z],-diff(htx,t)+diff(hxx,x)+diff(hxy,y)+diff(hxz,z))
#plt.plot(time,hp)
#temp1=diff(hyy,z,z)-2*diff(hyz,y,z)+diff(hzz,y,y)
#temp2=-diff(hty,t,y)-diff(htz,t,z)+diff(hxy,x,y)+diff(hxz,x,z)+diff(hyy,y,y)+2*diff(hyz,y,z)+diff(hzz,z,z)
#temp3=diff(htt,y,y)+diff(htt,z,z)-diff(hxx,y,y)-diff(hxx,z,z)-diff(hyy,y,y)-diff(hyy,z,z)-diff(hzz,y,y)-diff(hzz,z,z)
j1y=-diff(hyy,z)+diff(hyz,y)
#print(j1y)
j1z=-diff(hyz,z)+diff(hzz,y)
j2y=diff(htz,t)-diff(hxz,x)-diff(hyz,y)-diff(hzz,z)
j2z=-diff(hty,t)+diff(hxy,x)+diff(hyy,y)+diff(hyz,z)
j3y=diff(trace,z)/2
j3z=-diff(trace,y)/2
#Sapprox=7*(-diff(j1y,z)+diff(j1z,y))/6
#print(-diff(j1y,z)+diff(j1z,y))
S1=lambdify([t,x,y,z],(diff(j1y,z)-diff(j1z,y)))
print("one")
#S2=lambdify([t,x,y,z],diff(j2y,z)-diff(j2z,y))
#print(-diff(j2y,z)+diff(j2z,y))

#S3=lambdify([t,x,y,z],diff(j3y,z)-diff(j3z,y))
S=-diff(j1y+j2y+j3y,z)+diff(j1z+j2z+j3z,y)
#print(S)
Stota=lambdify([t,x,y,z],S)
print("two")
#print("three")
#plt.figure()
size=20
rarr=linspace(0,1,size)
tharr=linspace(0,2*pi,size)
rmesh,thmesh=np.meshgrid(rarr,tharr)
xmesh=rmesh*np.cos(thmesh)
ymesh=rmesh*np.sin(thmesh)
zp=[0.1,-0.5,1.75,-3]
#xp=0.4
#yp=0.3
k=50000


def Squad(zp):
    arr=[]
    arr2=[]
    for k in range(50000,200,-200):
        print(k)
        for i in range(size):
            for j in range(size):
                Smesh=np.zeros(2)
                Smesh[0]+=S1(c*time[k],xmesh[i,j],ymesh[i,j],zp)
                Smesh[1]+=Stota(c*time[k],xmesh[i,j],ymesh[i,j],zp)
                #Smesh[2]+=S3(c*time[k],xmesh[i,j],ymesh[i,j],zp)
        arr.append(np.sqrt(Smesh[1]**2/size**2))
        arr2.append(np.sqrt(Smesh[0]**2/size**2))
        

    return arr,arr2

#ani=FuncAnimation(fig,animate,range(50000,5000,-100),interval=10)
#plt.plot(S1(flip(c*time[50000:5000:-100]),xp,yp,zp)+S2(flip(c*time[50000:5000:-100]),xp,yp,zp)+S3(flip(c*time[50000:5000:-100]),xp,yp,zp))
#plt.plot(S1(flip(c*time[50000:5000:-100]),xp,yp,zp)/6)
fig, ax=plt.subplots(1,len(zp))
fig2, ax2=plt.subplots(1,len(zp))
for i in range(len(zp)):
    print(i)
    Stot,S1quad=Squad(zp[i])
    ax[i].plot(Stot,label=r'$S$')
    ax[i].plot(7/6*np.array(S1quad),label=r"$\frac{7}{6}\,S_1$")
    ax2[i].plot(np.array(Stot)/np.array(S1quad),label=r"$\frac{S}{S_1}$")
    ax2[i].plot(7/6*np.ones(len(Stot)),label=r"$\frac{7}{6}$")
    ax[i].set_ylim(0,1e-30)
    ax[i].set_xlabel(r'Time (arb. unit)')
    ax2[i].set_ylim(0,2)
    ax2[i].set_xlabel(r'Time (arb. unit)')
    #ax2[i].legend()
    ax[i].set_title(r"z={}".format(zp[i]))
    ax2[i].set_title(r"z={}".format(zp[i]))
ax[-1].legend(bbox_to_anchor =(1.45, 1.00))
ax[0].set_ylabel(r'Quadratic mean of the source term ($T^{-1}m^{-2}$)')
ax2[-1].legend(bbox_to_anchor =(1.35, 1.00))

plt.show()