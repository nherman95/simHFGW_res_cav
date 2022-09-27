from lal import MSUN_SI,PC_SI
from lalsimulation import SimInspiralChooseTDWaveform,TaylorT4,TaylorT2,SimInspiralChooseFDWaveform,EccentricFD ,TaylorF2,SpinTaylorT4Fourier,TaylorF2NLTides,TaylorR2F4
from math import log,ceil,log2
from numpy import sqrt,trapz,heaviside,gradient,pi,linspace,array,outer,fft,arange,zeros,empty,real,imag,complex128,sum,exp,sin,cos,dot
import numpy as np
from scipy.special import jv, yv
from scipy.optimize import newton_krylov,brentq
from scipy.interpolate import CubicSpline
from scipy.integrate import romb
from time import time
from scipy.constants import c,G
import numpy as np
from sys import argv

mu0 = 4*pi*10**(-7) #T #Tm/A ou kg m A-2s-2
Q=1e5
def J0(x):
    return jv(0,x)

def Rk_EMmax(x,A,Rmax):
    return A*jv(0,x*Rmax) + yv(0,x*Rmax)
def Rk_EMmin(x,A,Rmin):
    return A*jv(0,x*Rmin) + yv(0,x*Rmin)


def F(x,Rmax,Rmin):
    return yv(0,x*Rmax)*jv(0,x*Rmin)-yv(0,x*Rmin)*jv(0,x*Rmax)

def alphaAk_Rk_EM(kmax,rmax,rmin): #marche que pour kmax = 5  et rmin=0.1 et rmax>0.5 car fait à la main

    if rmax >= 1 :
        a=7/(4*rmax) #point de départ de recherche des zéros qui va varier (=7/(4*Rmax)) car 7 = première pseudo-période de J(x), celle ci décroit quand x augmente (il faut mieux overshoot la pseudo-période)
        i=0 #numéro des zéros trouver dans L
        ftol = 1e-6 #tolérance de convergence vers 0
        zero0 = float(newton_krylov(lambda x:F(x,rmax,rmin), a, f_tol=ftol)) #premier 0
        a = zero0
        Alpha = [zero0] # initialisation premier 0
        Alist = []


        while len(Alpha) < kmax :
            zero = float(newton_krylov(lambda x:F(x,rmax,rmin), a, f_tol=ftol))
            if zero-Alpha[i] >= 10*ftol :
                Alpha.append (zero)
                i+=1 # on va tester les zéros suivants avec le zéro qui vient d'être ajouté
            a+=7/(2*rmax) # la période la plus petite des 2 est celle minimum de recherche de zéro de la fonction F #ok avec 11/(4*Rmax) et 7/(2*Rmax)

    if 0.5 <= rmax < 1 :
        a=10/(4*rmax) # ok avec 10/(4*Rmax)
        i=0
        ftol = 1e-6
        zero0 = float(newton_krylov(lambda x:F(x,rmax,rmin), a, f_tol=ftol))
        a = zero0
        Alpha = [zero0]
        Alist = []



        while len(Alpha) < kmax :
            zero = float(newton_krylov(lambda x:F(x,rmax,rmin), a, f_tol=ftol))
            if zero-Alpha[i] >= 10*ftol :
                Alpha.append (zero)
                i+=1
            a+=(7.5/(2*rmax)) # ok avec (7.5/(2*Rmax))



    def G0(A):
        return Rk_EMmax(Alpha[0],A,rmax)
    Alist.append(brentq(G0,-1000,1000)) #super rapide (dicho sur droite)
    def G1(A):
        return Rk_EMmax(Alpha[1],A,rmax)
    Alist.append(brentq(G1,-1000,1000))
    def G2(A):
        return Rk_EMmax(Alpha[2],A,rmax)
    Alist.append(brentq(G2,-1000,1000))
    def G3(A):
        return Rk_EMmax(Alpha[3],A,rmax)
    Alist.append(brentq(G3,-1000,1000))
    def G4(A):
        return Rk_EMmax(Alpha[4],A,rmax)
    Alist.append(brentq(G4,-1000,1000))
    AlphaAlist =[]
    for j in range(0,len(Alpha)):
        AlphaAlist.append([Alpha[j],Alist[j]])
    return array(AlphaAlist) #liste de tuples (alpha_k, Ak)


def cav_parameters(cavity,K,l,Rmax,Rmin=0.1):
    if(cavity=='TM'):
        alpha_k=array([brentq(J0,i,i+1) for i in range(2,3*K,3)])/Rmax
        Ik=Rmax/alpha_k*jv(1,alpha_k*Rmax)
        normk=pi*0.5*(Rmax)**2*(jv(0,Rmax*alpha_k)**2 + jv(1,Rmax*alpha_k)**2)
    elif(cavity=='TEM'):
        alist=alphaAk_Rk_EM(K,Rmax,Rmin); alpha_k=array(alist[:,0])/l; Ak=array(alist[:,1])
        Ik=Ak/alpha_k*(Rmax*jv(1,alpha_k*Rmax)-Rmin*jv(1,alpha_k*Rmin)) + 1.0/alpha_k*(Rmax*yv(1,alpha_k*Rmax)-Rmin*yv(1,alpha_k*Rmin))
        normk=pi*0.5*(Rmax)**2*(Ak**2*(jv(0,Rmax*alpha_k)**2 + jv(1,Rmax*alpha_k)**2) + yv(0,Rmax*alpha_k)**2 + yv(1,Rmax*alpha_k)**2 + 2*Ak*(jv(0,Rmax*alpha_k)*yv(0,Rmax*alpha_k) + jv(1,Rmax*alpha_k)*yv(1,Rmax*alpha_k))) - pi*0.5*(Rmin)**2*(Ak**2*(jv(0,Rmin*alpha_k)**2 + jv(1,Rmin*alpha_k)**2) + yv(0,Rmin*alpha_k)**2 + yv(1,Rmin*alpha_k)**2 + 2*Ak*(jv(0,Rmin*alpha_k)*yv(0,Rmin*alpha_k) + jv(1,Rmin*alpha_k)*yv(1,Rmin*alpha_k)))
    return alpha_k,Ik,normk

def fct_bode(fspan,cavity,Rmax,K=5,l=1,B0=5):
    alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)
    testan=0
    for k in arange(len(alphak)):
        testan+=Ik[k]**2/normk[k]/sqrt((4*pi**2*fspan**2-alphak[k]**2*c**2/l**2)**2+(2*pi*fspan*alphak[k]*c/l)**2/Q**2)
    return c*2*sqrt(2)*pi**2*B0**2*l**2*4*pi**2*fspan**2/mu0*abs(sin(pi*fspan*l/c))*testan

def fchirp(f,M):
    A=sqrt(5/24)/pi**(2/3)
    chirpmass=M*MSUN_SI**(6/5)/(2*M*MSUN_SI)**(1/5)
    dist=1e9*PC_SI
    return A*c/dist*(G*chirpmass/c**3)**(5/6)/(2*pi*f)**(7/6)

def RMS_compute(func,fmin,fmax,cavity,Rmax,*args):
    # npts=1e7/211200000*(fmax-fmin)
    dpts=1e7
    # df=(fmax-fmin)*npts/dpts
    # num=round(npts/dpts)
    # #print(num)
    # res=0
    # for i in range(num):
        # print(i)
    ftest=linspace(fmin,fmax,2**round(log2(dpts))+1)
    a = (4*pi*fct_bode(ftest,cavity,Rmax)*func(ftest,*args))**2
    # plt.loglog(a)
    # plt.show()
    res=romb(a,ftest[1]-ftest[0])
    #print(res)
    return sqrt(res)

di=float(argv[1])
ds=float(argv[2])
num=int(argv[3])
Rmax=float(argv[4])
Rid=int(argv[5])
B0=5
l=1
cavity='TEM'
K=5
alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)
Mexp=np.loadtxt('mass.in')

Masse = 10**(-1*Mexp)

Deno = np.linspace(di,ds,num)
RMStemp=np.empty(len(Deno))
rep=np.empty(len(Mexp))
for j in range(len(Mexp)):
    fisco=2200/Masse[j]
    for i in range(len(Deno)):
        fmin=fisco/Deno[i]
        RMStemp[i]=RMS_compute(fchirp,fmin,fisco,cavity,Rmax,Masse[j])
    rep[j]=np.median(RMStemp)
title="results/simRay_{0:04d}.txt".format(Rid)
file = open(title, "w+")
np.savetxt(file,rep)
file.close()

