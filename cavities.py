from lal import MSUN_SI,PC_SI
from lalsimulation import SimInspiralChooseTDWaveform,TaylorT4,TaylorT2,SimInspiralChooseFDWaveform,EccentricFD ,TaylorF2,SpinTaylorT4Fourier,TaylorF2NLTides,TaylorR2F4
from math import log,ceil,log2
from numpy import sqrt,trapz,heaviside,gradient,pi,linspace,array,outer,fft,arange,zeros,empty,real,imag,complex128,sum,exp,sin,cos,dot
import numpy as np
from scipy.special import jv, yv,struve,jn_zeros
from scipy.optimize import root,fminbound
from scipy.interpolate import CubicSpline
from scipy.integrate import romb
from time import time
from scipy.constants import c,G
from sympy import besselj
from harmosc import harmosc_solver
from mpmath import meijerg,mpf,convert,hyp1f2


mu0 = 4*pi*10**(-7) #T #Tm/A ou kg m A-2s-2
Q=1e5
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
    zero0=root(lambda x: F(x,rmax,rmin),zeroj[0])
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
        zerok=fminbound(lambda x: F(x,rmax,rmin)**2,x1,x2)
        Alpha.append(zerok)
        Alist.append(-yv(1,zerok*rmin)/jv(1,zerok*rmin))
    AlphaAlist =[]
    for j in range(0,len(Alpha)):
        AlphaAlist.append((Alpha[j],Alist[j]))
    return array(AlphaAlist) #liste de tuples (alpha_k, Ak)


def cav_parameters(cavity,K,l,Rmax,Rmin=0.1):
    if(cavity=='TM'):
        alpha_k=jn_zeros(1,5)/Rmax
        Ik=-(jv(0,alpha_k*Rmax)-1)/alpha_k
        norm=1/sqrt(2/l/pi/alpha_k**2/jv(0,alpha_k*Rmax)**2)
    elif(cavity=='TEM'):
        alist=alphaAk_Rk_EM(K,Rmax,Rmin); alpha_k=array(alist[:,0])/l; Ak=array(alist[:,1])
        Ik=(Ak-Ak*jv(0,alpha_k*Rmax)-yv(0,alpha_k*Rmax))/alpha_k-(Ak-Ak*jv(0,alpha_k*Rmin)-yv(0,alpha_k*Rmin))/alpha_k
        norm=np.empty(5)
        for i in range(len(alpha_k)):
            tmp1=0.5*(-(Ak[i]**2+1)*hyp1f2(0.5,1,2,convert(str(-alpha_k[i]**2*Rmin**2)[1:-1]))-2/sqrt(pi)*meijerg([[],[-0.5,0.5]],[[-1,0,0],[-0.5]],convert(str(alpha_k[i]*Rmin)[1:-1]),0.5)-2*Ak[i]*jv(0,alpha_k[i]*Rmin)*yv(0,alpha_k[i]*Rmin)-2*Ak[i]*jv(1,alpha_k[i]*Rmin)*yv(1,alpha_k[i]*Rmin)+Ak[i]**2+1)
            #tmp1=1/2*Rmin**2*(2*Ak[i]/sqrt(pi)*meijerg([[0,0.5],[-0.5]],[[0,1],[-1,-1,-0.5]],convert(str(alpha_k[i]*Rmin)[1:-1]),0.5)+Ak[i]**2*(jv(1,alpha_k[i]*Rmin)**2-jv(0,alpha_k[i]*Rmin)*jv(2,alpha_k[i]*Rmin))+yv(1,alpha_k[i]*Rmin)**2-yv(0,alpha_k[i]*Rmin)*yv(2,alpha_k[i]*Rmin))
            #tmp2=1/2*Rmax**2*(2*Ak[i]/sqrt(pi)*meijerg([[0,0.5],[-0.5]],[[0,1],[-1,-1,-0.5]],convert(str(alpha_k[i]*Rmax)[1:-1]),0.5)+Ak[i]**2*(jv(1,alpha_k[i]*Rmax)**2-jv(0,alpha_k[i]*Rmax)*jv(2,alpha_k[i]*Rmax))+yv(1,alpha_k[i]*Rmax)**2-yv(0,alpha_k[i]*Rmax)*yv(2,alpha_k[i]*Rmax))
            tmp2=0.5*(-(Ak[i]**2+1)*hyp1f2(0.5,1,2,convert(str(-alpha_k[i]**2*Rmax**2)[1:-1]))-2/sqrt(pi)*meijerg([[],[-0.5,0.5]],[[-1,0,0],[-0.5]],convert(str(alpha_k[i]*Rmax)[1:-1]),0.5)-2*Ak[i]*jv(0,alpha_k[i]*Rmax)*yv(0,alpha_k[i]*Rmax)-2*Ak[i]*jv(1,alpha_k[i]*Rmax)*yv(1,alpha_k[i]*Rmax)+Ak[i]**2+1)
            norm[i]=sqrt((tmp2-tmp1)*l/2*pi)
        #norm=1/sqrt(2/l/pi/alpha_k**2/(Ak**2*(jv(0,alpha_k*Rmax)**2)+(yv(0,alpha_k*Rmax)**2)))#-(Ak*jv(0,alpha_k*Rmin)+yv(0,alpha_k*Rmin))**2))
    return alpha_k.flatten(),Ik.flatten(),norm.flatten()

def gw_binary_time_signal(M,dist,f_min,deltaT):
    distance = dist * PC_SI 
    m = M * MSUN_SI; other=0.0; params=None
    hp,hx=SimInspiralChooseTDWaveform( m,  m,  other,  other,  other, other,  other,  other,  distance,  other,  other,  other, other,  other,  deltaT,  f_min, other, params, TaylorT4)
    tmax=-float(hp.epoch); hp=hp.data.data;t=linspace(0.0,tmax,len(hp));hx=hx.data.data
    return hp,hx,t

def gw_binary_freq_signal(M,dist,f_min,f_max,deltaF):
    distance = dist * PC_SI
    m = M * MSUN_SI; other=0.0; params=None
    TFhp,_=SimInspiralChooseFDWaveform( m,  m,  other,  other,  other, other,  other,  other,  distance,  other,  other,  other, other,  other,  deltaF,  f_min, f_max, other, params, TaylorF2 )
    TFhp=2*deltaF*TFhp.data.data
    omega = 2*pi*linspace(0,f_max,len(TFhp))
    return TFhp,omega



def time_approach(hp,t,cavity,Rmax,K=5,l=1,B0=5):
    t0 = time()
    alphak,Ik,norm=cav_parameters(cavity,K,l,Rmax)
    T=c/l*t ; dt=t[1]-t[0] ; dT=T[1]-T[0] ; tmax=t[-1]
    hGW=max(abs(hp)) ; hp=hp/hGW
    hinter=CubicSpline(T,hp); nbrzeros=200
    TD=T-1; TD[TD<0]=0 ; scoef=pi*Ik/norm 
    sint=hinter(T,1)-hinter(TD,1)
    sk=outer(scoef,sint)
    freq=2.*pi*fft.rfftfreq(len(T)+nbrzeros,dT); coef=[]
    for k in arange(len(alphak)):
       coef.append((alphak[k]**2-freq**2+1j*freq*alphak[k]/Q))
    f = c/l*fft.rfftfreq(len(sk[0])+200,dT)
    transfo=2./(len(T)+nbrzeros)*fft.rfft(sk,axis=1,n=len(T)+nbrzeros)/coef

    # rescaling sampling if necessary
    if (1.0/(2.0*dt)<alphak[-1]*c/l/2.0/pi):
        fs=3*alphak[-1]/2.0/pi*c/l
        nvec=2**(int(ceil(log(fs*t[-1],2))))
        T=linspace(0,t[-1],nvec)*c/l
        dT=T[1]-T[0]
        print(len(T))

    harmosc_solver.solver(T,alphak,freq,real(transfo),imag(transfo),Q)
    
    dE_t=(-7/3*pi*B0**2/mu0*hGW)*dot(Ik/norm,harmosc_solver.sol)
    harmosc_solver.sol=None
    
    tdE_t=linspace(0,tmax,len(dE_t))
    dP_t = gradient(dE_t,tdE_t[1]-tdE_t[0],edge_order=2)
    rmsP_t = sqrt(sum(abs(dP_t)**2)/len(dP_t))#sqrt(trapz(dP_t**2,tdE_t)/tmax)
    print('time approach done in (s)')
    print(time()-t0)
    return rmsP_t,dP_t


def freq_appoach(TFhp,omega,cavity,Rmax,recomp=0,K=5,l=1,B0=5):
    alphak,Ik,norm=cav_parameters(cavity,K,l,Rmax)
    deltaF=1/2/pi*(omega[1]-omega[0])
    if (omega[-1]<alphak[-1]*c):
        omega=np.arange(omega[0],alphak[-1]*c+100*pi*deltaF,omega[1]-omega[0])
        TFhp=np.pad(TFhp,(0,len(omega)-len(TFhp)),constant_values=(0,0))
    TFsk = 7/3*pi*B0*c*outer(Ik/norm,omega*TFhp*sin(omega*l/2/c))
    coef=[]
    for k in arange(len(alphak)):
        coef.append(((alphak[k]*(c/l))**2-omega**2+1j*omega*alphak[k]/Q*c))
    Ak = array(-sum(real(TFsk/coef),axis=1))
    Bk = array((-alphak/2/Q*c*Ak+sum(omega*imag(TFsk/coef),axis=1))/(alphak*c/l)/sqrt(1-1/2/Q))
    dirack=[]
    for k in range (0,len(alphak)):
        a=Ak[k]/(-alphak[k]/2/Q*c*l+1j*(omega-alphak[k]*c/l*sqrt(1-1/2/Q)))
        b=1j*Bk[k]/(-alphak[k]/2/Q*c*l+1j*(omega-alphak[k]*c/l*sqrt(1-1/2/Q)))
        dirack.append(a+b)
    TFbk = TFsk/coef + array(dirack)
    TFdE_f_p = 2*pi*B0/mu0*(dot(Ik/norm,TFbk))#*heaviside(-omega/2/pi+fisco,0.5)*heaviside(omega/2/pi-f_min,0.5)
    TFdP_f_p = 1j*omega*TFdE_f_p
    rmsP_f = sqrt(sum(abs(TFdP_f_p)**2)/len(TFdP_f_p))
    dP_f = fft.irfft(2j*omega*pi*B0/mu0*(dot(Ik/norm,TFsk/coef)))*len(TFdP_f_p)
    tdE_f = linspace(0,1/deltaF,len(dP_f))
    if (recomp==1):
        freq=alphak*c*sqrt(1-1/2/Q)
        for i in range(len(tdE_f)):
            dP_f[i]=dP_f[i]+(dot((Ak)*exp(-alphak/2/Q*c*tdE_f[i]),-freq*sin((freq-deltaF)*tdE_f[i]))+dot((Bk)*exp(-alphak/2/Q*c*tdE_f[i]),freq*cos((freq-deltaF)*tdE_f[i])))*pi*B0/mu0-(dot(alphak/2/Q*c*exp(-alphak/2/Q*c*tdE_f[i])*(Ak),cos((freq-deltaF)*tdE_f[i]))+dot(alphak/2/Q*c*exp(-alphak/2/Q*c*tdE_f[i])*(Bk),sin((freq-deltaF)*tdE_f[i])))*2*pi*B0/mu0
    return rmsP_f,TFdP_f_p,omega,dP_f,tdE_f

def fct_bode(fspan,cavity,Rmax,K=5,l=1,B0=5):
    alphak,Ik,norm=cav_parameters(cavity,K,l,Rmax)
    testan=0
    for k in arange(len(alphak)):
        testan+=2*Ik[k]**2/norm[k]**2/sqrt((4*pi**2*fspan**2-alphak[k]**2*c**2/l**2)**2 +(alphak[k]*2*pi*fspan*c/2/Q)**2)
    return c*7*pi**2*sqrt(2)*B0**2*4*pi**2*fspan**2/6/mu0*abs(sin(pi*fspan*l/c))*testan

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
    a = 2*pi*(fct_bode(ftest,cavity,Rmax)*func(ftest,*args))**2
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

def sigmoid(x,fcut):
    return 0.5+0.5*np.tanh((fcut-x)/fcut*20)

def Sh_toy(f,omegaGW,fcut):
    omegaGW=1e-10
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*heaviside(-f+fcut,0.5)
    return sqrt(Sf/2)
def Sh_toyfm1(f,omegaGW,fcut):
    omegaGW=1e-10*1e7/f
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*heaviside(-f+fcut,0.5)
    return sqrt(Sf/2)
def Sh_toyff1(f,omegaGW,fcut):
    omegaGW=1e-10/1e7*f
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*heaviside(-f+fcut,0.5)
    return sqrt(Sf/2)
def Sh_toyfm2(f,omegaGW,fcut):
    omegaGW=1e-10*1e7/f**2
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*heaviside(-f+fcut,0.5)
    return sqrt(Sf/2)
def Sh_toyff2(f,omegaGW,fcut):
    omegaGW=1e-10/1e7*f**2
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*heaviside(-f+fcut,0.5)
    return sqrt(Sf/2)

def SSh_toy(f,omegaGW,fcut):
    omegaGW=1e-10
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*sigmoid(f,fcut)
    return sqrt(Sf/2)
def SSh_toyfm1(f,omegaGW,fcut):
    omegaGW=1e-10*1e7/f
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*sigmoid(f,fcut)
    return sqrt(Sf/2)
def SSh_toyff1(f,omegaGW,fcut):
    omegaGW=1e-10/1e7*f
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*sigmoid(f,fcut)
    return sqrt(Sf/2)
def SSh_toyfm2(f,omegaGW,fcut):
    omegaGW=1e-10*1e7/f**2
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*sigmoid(f,fcut)
    return sqrt(Sf/2)
def SSh_toyff2(f,omegaGW,fcut):
    omegaGW=1e-10/1e7*f**2
    H0=70/1e3/PC_SI
    Sf=3*H0**2/4/pi**2/f**3*omegaGW*sigmoid(f,fcut)
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

   
