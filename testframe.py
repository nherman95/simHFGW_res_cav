from cavities import *
from sympy.abc import t,x,y,z
from sympy import diff,init_printing,simplify
from sympy.utilities.lambdify import lambdify
import sympy
import numpy as np
from numpy import linspace,flip
from math import log,factorial
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

mu=1e-10/2e-5*MSUN_SI*G/c**2
m=2e-5*MSUN_SI*G/c**2
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
hpan2s=diff(hpan2s,t,t)
hxan2s=diff(hxan2s,t,t)
print(hpan2s)
Pp=1/6*z**2*hpan2s
Px=1/6*z**2*hxan2s
Qp=1/3*z**2*hpan2s
Qx=1/3*z**2*hxan2s
for n in range (3,8):
    print(n)
    hpan2s=diff(hpan2s,t)
    Pp=Pp+(n-1)/factorial(n+1)*z**n*hpan2s
    Px=Px+(n-1)/factorial(n+1)*z**n*hxan2s
    Qp=Qp+(n)/factorial(n+1)*z**n*hpan2s
    Qx=Qx+(n)/factorial(n+1)*z**n*hxan2s

hxx=Pp
hyy=-Pp
hxy=Px
hxz=-1/z*(x*Pp+y*Px)
hyz=-1/z*(x*Px-y*Pp)
hzz=1/z**2*((x**2-y**2)*Pp+2*x*y*Px)
htt=1/z**2*((x**2-y**2)*(2*Qp-Pp)+2*x*y*(2*Qx-Px))
htx=-1/z*(x*Qp+y*Qx)
hty=-1/z*(x*Qx-y*Qp)
htz=1/z**2*((x**2-y**2)*(Qp)+2*x*y*(Qx))
trace=-htt+hxx+hyy+hzz
#dondert=lambdify([t,x,y,z],-diff(htt,t)+diff(htx,x)+diff(hty,y)+diff(htz,z))
#donderx=lambdify([t,x,y,z],-diff(htx,t)+diff(hxx,x)+diff(hxy,y)+diff(hxz,z))
#plt.plot(time,hp)
temp1=diff(hyy,z,z)-2*diff(hyz,y,z)+diff(hzz,y,y)
temp2=-diff(hty,t,y)-diff(htz,t,z)+diff(hxy,x,y)+diff(hxz,x,z)+diff(hyy,y,y)+2*diff(hyz,y,z)+diff(hzz,z,z)
temp3=diff(htt,y,y)+diff(htt,z,z)-diff(hxx,y,y)-diff(hxx,z,z)-diff(hyy,y,y)-diff(hyy,z,z)-diff(hzz,y,y)-diff(hzz,z,z)
j1y=-diff(hyy,z)+diff(hyz,y)
j1z=-diff(hyz,z)+diff(hzz,y)
j2y=diff(htz,t)-diff(hxz,x)-diff(hyz,y)-diff(hzz,z)
j2z=diff(hty,t)+diff(hxy,x)+diff(hyy,y)+diff(hyz,z)
j3y=diff(trace,z)/2
j3z=-diff(trace,y)/2
Sapprox=-diff(Pp,z,z)-diff(Pp,z)/z-diff(Qp,t)/z-3*Pp/z**2+2*Qp/z**2
S1=lambdify([t,x,y,z],(diff(j1y,z)-diff(j1z,y)))
print("one")
S2=lambdify([t,x,y,z],diff(j2y,z)-diff(j2z,y))
print("two")
S3=lambdify([t,x,y,z],diff(j3y,z)-diff(j3z,y))
Stota=lambdify([t,z],-Sapprox)
print("three")
#plt.figure()
size=10
rarr=linspace(0,1,size)
tharr=linspace(0,2*pi,size)
rmesh,thmesh=np.meshgrid(rarr,tharr)
xmesh=rmesh*np.cos(thmesh)
ymesh=rmesh*np.sin(thmesh)
zp=-0.5 
xp=0.4
yp=0.3
Smesh=np.empty((4,size,size))
fig, ax = plt.subplots(1,3)
k=50000
for i in range(size):
        for j in range(size):
            Smesh[0,i,j]=S1(c*time[k],xmesh[i,j],ymesh[i,j],zp)
            Smesh[1,i,j]=S2(c*time[k],xmesh[i,j],ymesh[i,j],zp)
            Smesh[2,i,j]=S3(c*time[k],xmesh[i,j],ymesh[i,j],zp)
            Smesh[3,i,j]=Stota(c*time[k],zp)
mesh_1=ax[0].pcolormesh(xmesh, ymesh, np.abs(Smesh[0,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
mesh_2=ax[1].pcolormesh(xmesh, ymesh, np.abs(Smesh[3,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
#mesh_3=ax[2].pcolormesh(xmesh, ymesh, np.abs(Smesh[2,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
mesh_4=ax[2].pcolormesh(xmesh, ymesh, np.abs(Smesh[2,:,:]+Smesh[0,:,:]+Smesh[1,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
ax[0].set_box_aspect(1)
ax[1].set_box_aspect(1)
ax[2].set_box_aspect(1)
#ax[3].set_box_aspect(1)
fig.colorbar(mesh_2, ax=ax)


def animate(k):
    for i in range(size):
        for j in range(size):
            Smesh[0,i,j]=S1(c*time[k],xmesh[i,j],ymesh[i,j],zp)
            Smesh[1,i,j]=S2(c*time[k],xmesh[i,j],ymesh[i,j],zp)
            Smesh[2,i,j]=S3(c*time[k],xmesh[i,j],ymesh[i,j],zp)
            Smesh[3,i,j]=Stota(c*time[k],zp)
    mesh_1=ax[0].pcolormesh(xmesh, ymesh, np.abs(7*Smesh[0,:,:]/6),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
    mesh_2=ax[1].pcolormesh(xmesh, ymesh, np.abs(Smesh[3,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
    #mesh_3=ax[2].pcolormesh(xmesh, ymesh, np.abs(Smesh[2,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
    mesh_4=ax[2].pcolormesh(xmesh, ymesh, np.abs(Smesh[2,:,:]+Smesh[0,:,:]+Smesh[1,:,:]),cmap='coolwarm',norm=colors.LogNorm(vmin=1e-33,vmax=1e-30))
    ax[0].set_box_aspect(1)
    ax[1].set_box_aspect(1)
    ax[2].set_box_aspect(1)
    #ax[3].set_box_aspect(1)
    # fig.colorbar(mesh_1, ax=ax[0])
    # fig.colorbar(mesh_3, ax=ax[2])
    return 

ani=FuncAnimation(fig,animate,range(50000,5000,-100),interval=10)
#plt.plot(S1(flip(c*time[50000:5000:-100]),xp,yp,zp)+S2(flip(c*time[50000:5000:-100]),xp,yp,zp)+S3(flip(c*time[50000:5000:-100]),xp,yp,zp))
#plt.plot(S1(flip(c*time[50000:5000:-100]),xp,yp,zp)/6)
#plt.plot(Stota(flip(c*time[50000:5000:-100]),zp))
plt.show()