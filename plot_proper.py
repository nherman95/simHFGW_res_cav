import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scis
import scipy.optimize as scio
import scipy.integrate as scii
import scipy.interpolate as sciinter
from time import time
from matplotlib import rcParams
from scipy.special import jv, yv,struve,jn_zeros
from numpy import sqrt,trapz,heaviside,gradient,pi,linspace,array,outer,fft,arange,zeros,empty,real,imag,complex128,sum,exp,sin,cos,dot
from scipy.optimize import newton_krylov,brentq


L = 1 #m
#Rmax = 5 #m
Rmin = 0.1 #m
B0 = 5 #T
mu0 = 4*np.pi*10**(-7) #Tm/A ou kg m A-2s-2
c = 299792458.0 #m/s #m/s


## Définitions des paramètres de la cavité (alphak, Rrondk Ik, pour TM et TEM) 2/5
# pour les conditions limites(r fixé à Rmax ou Rmin)
def J1(x):
    return jv(1,x)

def Rk_EMmax(x,A,Rmax):
    return A*jv(1,x*Rmax) + yv(1,x*Rmax)
def Rk_EMmin(x,A,Rmin):
    return A*jv(1,x*Rmin) + yv(1,x*Rmin)


def F(x,Rmax,Rmin):
    return yv(1,x*Rmax)*jv(1,x*Rmin)-yv(1,x*Rmin)*jv(1,x*Rmax)

def alphaAk_Rk_EM(kmax,rmax,rmin): #marche que pour kmax = 5  et rmin=0.1 et rmax>0.5 car fait à la main
    zeroj=jn_zeros(1,kmax)/rmax
    zero0=scio.root(lambda x: F(x,rmax,rmin),zeroj[0])
    Alpha=[zero0.x]
    Alist=[-yv(1,zero0.x*rmin)/jv(1,zero0.x*rmin)]
    k=0
    for k in range(1,kmax):
        if(rmax<0.15):
            x1=Alpha[k-1]*1.001
            x2=x1+(zeroj[k]-zeroj[k-1])*3
            # while(F(x1,rmax,rmin)*F(x2,rmax,rmin)>0):
            #     x2=x1+0.1
        elif(0.15<=rmax<0.2):
            x1=Alpha[k-1]*1.01
            x2=x1+(zeroj[k]-zeroj[k-1])*3
            # while(F(x1,rmax,rmin)*F(x2,rmax,rmin)>0):
            #     x2=x1+0.1
        else:
            x1=Alpha[k-1]*1.1
            x2=x1+(zeroj[k]-zeroj[k-1])*2
        zerok=scio.fminbound(lambda x: F(x,rmax,rmin)**2,x1,x2)
        Alpha.append(zerok)
        Alist.append(-yv(1,zerok*rmin)/jv(1,zerok*rmin))
    AlphaAlist =[]
    for j in range(0,len(Alpha)):
        AlphaAlist.append([Alpha[j],Alist[j]])
    return array(AlphaAlist) #liste de tuples (alpha_k, Ak)



t0 = time()
Rmax = 0.15
R =[]
F0_EM=[]
F1_EM=[]
F2_EM=[]
F3_EM=[]
F4_EM=[]
F0_M=[]
F1_M=[]
F2_M=[]
F3_M=[]
F4_M=[]

rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 10, 10
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 36
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 8
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
rcParams['text.usetex'] = True

while(Rmax < 10):
    print(Rmax)
    AlphaA_EM = alphaAk_Rk_EM(5,Rmax,Rmin)

    Alpha_M = jn_zeros(1,5)/Rmax

    R.append(Rmax)

    F0_EM.append(AlphaA_EM[0][0]*3*10**8/(2*np.pi))
    F1_EM.append(AlphaA_EM[1][0]*3*10**8/(2*np.pi))
    F2_EM.append(AlphaA_EM[2][0]*3*10**8/(2*np.pi))
    F3_EM.append(AlphaA_EM[3][0]*3*10**8/(2*np.pi))
    F4_EM.append(AlphaA_EM[4][0]*3*10**8/(2*np.pi))

    F0_M.append(Alpha_M[0]*3*10**8/(2*np.pi))
    F1_M.append(Alpha_M[1]*3*10**8/(2*np.pi))
    F2_M.append(Alpha_M[2]*3*10**8/(2*np.pi))
    F3_M.append(Alpha_M[3]*3*10**8/(2*np.pi))
    F4_M.append(Alpha_M[4]*3*10**8/(2*np.pi))

    if 0.13 > Rmax:
        Rmax += 0.005
    if 0.2 > Rmax >= 0.13:
        Rmax += 0.005
    if 0.5 > Rmax >= 0.2:
        Rmax += 0.01
    if 1 > Rmax >= 0.5:
        Rmax += 0.02
    if Rmax >= 1:
        Rmax += 0.1

Rmax = 5 # on réinitialise Rmax
t1 = time()
print(t1-t0)
plt.figure()
ax=plt.gca()
plt.loglog(R,F0_EM,'--',color="b",linewidth=1)
plt.loglog(R,F1_EM,'--',color="r",linewidth=1)
plt.loglog(R,F2_EM,'--',color="g",linewidth=1)
plt.loglog(R,F3_EM,'--',color="k",linewidth=1)
plt.loglog(R,F4_EM,'--',color="m",linewidth=1)

plt.loglog(R,F0_M,color="b",linewidth=1,label="k = 0")
plt.loglog(R,F1_M,color="r",linewidth=1,label="k = 1")
plt.loglog(R,F2_M,color="g",linewidth=1,label="k = 2")
plt.loglog(R,F3_M,color="k",linewidth=1,label="k = 3")
plt.loglog(R,F4_M,color="m",linewidth=1,label="k = 4")
plt.xlim(0.1,10)
ax.set_box_aspect(1)
plt.legend()
plt.xlabel(r"Outer radius (m)")
plt.ylabel(r"Frequency (Hz)")
plt.grid(True, which="both")
plt.savefig("proper.png")
