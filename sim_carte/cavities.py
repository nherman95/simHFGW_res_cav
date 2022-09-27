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
from multiprocess import Pool
import matplotlib.pyplot as plt
from harmosc import harmosc_solver


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

def gw_binary_time_signal(M,dist,f_min,deltaT):
    distance = dist * PC_SI 
    m = M * MSUN_SI; other=0.0; params=None
    hp,_=SimInspiralChooseTDWaveform( m,  m,  other,  other,  other, other,  other,  other,  distance,  other,  other,  other, other,  other,  deltaT,  f_min, other, params, TaylorT4)
    tmax=-float(hp.epoch); hp=hp.data.data;t=linspace(0.0,tmax,len(hp))
    return hp,t

def gw_binary_freq_signal(M,dist,f_min,f_max,deltaF):
    distance = dist * PC_SI
    m = M * MSUN_SI; other=0.0; params=None
    TFhp,_=SimInspiralChooseFDWaveform( m,  m,  other,  other,  other, other,  other,  other,  distance,  other,  other,  other, other,  other,  deltaF,  f_min, f_max, other, params, EccentricFD )
    TFhp=2*deltaF*TFhp.data.data
    omega = 2*pi*linspace(0,f_max,len(TFhp))
    return TFhp,omega


def freq_appoach(TFhp,omega,fmin,fmax,cavity,Rmax,K=5,l=1,B0=5):

    alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)
    deltaF=1/2/pi*(omega[1]-omega[0])
    TFsk = -2*pi*B0*l**2*c*outer(Ik/normk,omega*TFhp*sin(omega*l/2/c))
    coef=[]
    for k in arange(len(alphak)):
        coef.append(((alphak[k]*(c/l))**2-omega**2+1j*omega/Q/c**2))
    Ak = array(-sum(real(TFsk/coef),axis=1))
    Bk = array((-alphak/2/Q/c*l*Ak+sum(omega*imag(TFsk/coef),axis=1))/(alphak*c/l)/sqrt(1+1/4/Q))
    dirack=[]
    for k in range (0,len(alphak)):
        a=Ak[k]/(-alphak[k]/2/Q/c*l+1j*(omega-alphak[k]*c/l))+1j*Bk[k]/(-alphak[k]/2/Q/c*l+1j*(omega-alphak[k]*c/l))
        dirack.append(a)
    TFbk = TFsk/coef + array(dirack)
    TFdE_f_p = 4*pi*B0/mu0*dot(Ik,TFbk)#*heaviside(-omega/2/pi+fisco,0.5)*heaviside(omega/2/pi-f_min,0.5)
    # plt.loglog(omega,abs(TFdE_f_p))
    # plt.show()
    TFdP_f_p = 1j*omega*TFdE_f_p

    dP_f = fft.irfft(TFdP_f_p)*len(TFdP_f_p)
    tdE_f = linspace(0,1/deltaF,len(dP_f))

    rmsP_f = sqrt(trapz(abs(TFdP_f_p)**2,omega)/omega[-1])
    return rmsP_f,TFdP_f_p,dP_f,tdE_f

def time_approach(hp,t,cavity,Rmax,K=5,l=1,B0=5):
    t0 = time()
    alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)
    T=c/l*t ; dt=t[1]-t[0] ; dT=T[1]-T[0] ; tmax=t[-1]
    hGW=max(abs(hp)) ; hp=hp/hGW
    hinter=CubicSpline(T,hp); nbrzeros=200
    TD=T-1; TD[TD<0]=0 ; scoef=pi*Ik/normk
    sint=hinter(T,1)-hinter(TD,1)
    sk=outer(scoef,sint)
    f_min=2200/25/1e-5/c
    freq=2.*pi*fft.rfftfreq(len(T)+nbrzeros,dT); coef=[]
    for k in arange(len(alphak)):
       coef.append((alphak[k]**2-freq**2+1j*freq*alphak[k]/Q))
    f = c/l*fft.rfftfreq(len(sk[0])+200,dT)
    transfo=2./(len(T)+nbrzeros)*fft.rfft(sk,axis=1,n=len(T)+nbrzeros)/coef*heaviside(f-f_min,0.5)

    # rescaling sampling if necessary
    if (1.0/(2.0*dt)<alphak[-1]*c/l/2.0/pi):
        fs=3*alphak[-1]/2.0/pi*c/l
        nvec=2**(int(ceil(log(fs*t[-1],2))))
        T=linspace(0,t[-1],nvec)*c/l
        dT=T[1]-T[0]
        print(len(T))

    harmosc_solver.solver(T,alphak,freq,real(transfo),imag(transfo),Q)
    
    print(harmosc_solver.sol[:,0])
    dE_t=(4*pi*B0**2/mu0*hGW)*dot(Ik,harmosc_solver.sol)
    print(dE_t[0])
    harmosc_solver.sol=None
    tdE_t=linspace(0,tmax,len(dE_t))
    dP_t = gradient(dE_t,tdE_t[1]-tdE_t[0],edge_order=2)
    print(dP_t[0])
    rmsP_t = sqrt(trapz(dP_t**2,tdE_t)/tmax)
    print('time approach done in (s)')
    print(time()-t0)
    return rmsP_t,dP_t

def fct_bode(fspan,cavity,Rmax,K=5,l=1,B0=5):
    alphak,Ik,normk=cav_parameters(cavity,K,l,Rmax)
    testan=0
    for k in arange(len(alphak)):
        testan+=Ik[k]**2/normk[k]/sqrt((4*pi**2*fspan**2-alphak[k]**2*c**2/l**2)**2+(2*pi*fspan*alphak[k]*c/l)**2/Q**2)
    return c*2*sqrt(2)*pi**2*B0**2*l**2*4*pi**2*fspan**2/mu0*abs(sin(pi*fspan*l/c))*testan

def f_RMS_mass_deno(ftest,M,cavity,Rmax):
    A=sqrt(5/24)/pi**(2/3)
    chirpmass=M*MSUN_SI**(6/5)/(2*M*MSUN_SI)**(1/5)
    dist=1e9*PC_SI
    rep=A*c/dist*(G*chirpmass/c**3)**(5/6)/(2*pi*ftest)**(7/6)
    a = (4*pi*fct_bode(ftest,cavity,Rmax)*rep)**2
    return a

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






## SGWB en fct de rmax
alpha_infl = 2#4/3 # modèle de ps qu'on veut en c*f^alpha (alpha compris entre 2 et -2 pour retour simu temporel en freq osef)
alpha1_pt = 2
alpha2_pt = 4
alpha_ph = 2 #2.4

fcut1 = 10**8 #Hz
fcut2 = 8*10**8
fcut3 = 10**9
fcut4 = 10**7
fcut5=8*10**12

hc_infl = 0.2*10**-31 #au niveau du cut
hc_pt = 0.1*10**-30
hc1_ph = 0.3*10**-31 #borne sup de la box
hc2_ph = 10**-33

def Sh_infl(f):
    return sqrt((hc_infl**2)*(f/fcut1)**(-alpha_infl)/f/2)*heaviside(-f+fcut1,0.5)

# def Sh1_pt(f):
#     return (hc_pt**2)*(f/fcut2)**(-alpha1_pt)/f
# def Sh2_pt(f):
#     return (hc_pt**2)*(f/fcut2)**(-alpha2_pt)/f

def Sh_pt(f):
    rep=np.empty(len(f))
    rep[f <= fcut2]=(hc_pt**2)*(f[f <= fcut2]/fcut2)**(-alpha1_pt)/f[f <= fcut2]
    rep[f>fcut2]=(hc_pt**2)*(f[f>fcut2]/fcut2)**(-alpha2_pt)/f[f>fcut2]
    rep[f>fcut5]=0
    return sqrt(rep/2)


# def Sh1_ph(f):
#     return (hc1_ph**2)*(f/fcut1)**(-alpha_ph)
# def Sh2_ph(f):
#     return (hc2_ph**2)*(f/fcut3)**(-alpha_ph)

def Sh_ph(f):
    rep=np.zeros(len(f))
    rep[f <= fcut3]=(hc2_ph**2)*(f[f<= fcut3]/fcut3)**(-alpha_ph)/f[f<= fcut3]
    rep[f <= fcut1]=(hc1_ph**2)*(f[ f <= fcut1]/fcut1)**(-alpha_ph)/f[ f <= fcut1]
    rep[f<= fcut4]=0
    return sqrt(rep/2)

def Sh_toy(f):
    omegaGW=1e-10
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*heaviside(-f+fcut1,0.5)
    return sqrt(Sf/2)

def get_sgwbs(freq):
    TFhp_infl = Sh_infl(freq)
    TFhp_pt=Sh_pt(freq)
    TFhp_ph=Sh_ph(freq)
    # TFhp_pt = []
    # for f in freq :
    #     if f <= fcut2:
    #         TFhp_pt.append(sqrt(Sh1_pt(f)))
    #     else :
    #         TFhp_pt.append(sqrt(Sh2_pt(f)))
    # TFhp_pt = array(TFhp_pt)
    # TFhp_ph = []
    # for f in freq :
    #     if f <= fcut4 :
    #         TFhp_ph.append(0)
    #     if fcut4 < f <= fcut1:
    #         TFhp_ph.append(sqrt(Sh1_ph(f)))
    #     if fcut1 < f <= fcut3:
    #         TFhp_ph.append(sqrt(Sh2_ph(f)))
    #     if f > fcut3:
    #         TFhp_ph.append(0)
    # TFhp_ph = array(TFhp_ph)
    return TFhp_infl,TFhp_ph,TFhp_pt

   
