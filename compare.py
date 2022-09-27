from scipy.optimize.optimize import fmin
from cavities import *
import matplotlib.pyplot as plt
cavity='TM'
Mexp=5
den_f=15
M=10**(-Mexp); #Msun
print(M)
fisco=2200/(M); #hz
f_min = fisco/den_f;# start frequency of inspiral
f_max=fisco

if -log(M) < 5:
    f_max = 10**(log(M)+5)*fisco


l = 1.0  ;rmax = 5.0 ; rmin = 0.1 #m
Rmax=rmax/l ; Rmin=rmin/l
B0 = 5.0 ; 
K=5; # number of radial modes

dist = 1e9; deltaT = 1.0/2**(int(round(log(fisco,2)))+2);
hp,hx,t=gw_binary_time_signal(M,dist,f_min,deltaT)
tmax=t[-1] ; dt=t[1]-t[0]
deltaF = 1/tmax
TFhp,omega=gw_binary_freq_signal(M,dist,f_min,f_max,deltaF)
# plt.loglog(omega,abs(TFhp))
# plt.loglog(2*pi*linspace(f_min,f_max),fchirp(linspace(f_min,f_max),M)*deltaF*2)
# plt.show()
if -log(M) < 5:
    TFhp = TFhp*heaviside(-omega/2/pi+fisco,0.5)

#TFhp1=interp1d(omega,TFhp,kind='nearest')
#omega1=2*pi*linspace(0,f_max,10*len(TFhp))
rmsP_t,dP_t=time_approach(hp,t,cavity,Rmax)
rmsP_f,TFdP_f,dP_f,t_f=freq_appoach(TFhp,omega,cavity,Rmax)
t=linspace(0,t[-1],len(dP_t))
print(dP_f[0])
#

#TFdE_t=fft.rfft(dE_t)*2/len(dE_t)
#freqdE_t = 2*pi*fft.rfftfreq(len(T),dT)
TFdP_t = fft.rfft(dP_t,n=len(dP_t))*2/len(dP_t)
freqdP_t = fft.rfftfreq(len(dP_t),t[1]-t[0])

# dP_t_fit = CubicSpline(t,dP_t)
# dP_f_fit = CubicSpline(t_f-0*dt,dP_f)
# diff = ((dP_t_fit(t)) - (dP_f_fit(t)))
# ecart = sqrt(trapz(diff[0:len(diff)-100]**2,t[0:len(diff)-100])/tmax)/((rmsP_f+rmsP_t)/2)
#print(ecart)
print(rmsP_t)
print(rmsP_f)
#print(sqrt(trapz(abs(TFdP_t)**2,2*pi*freqdP_t)/2/pi/freqdP_t[-1]))
#print(RMS_mass_deno(1e-5,25,cavity,Rmax))
# plt.loglog(omega/2/pi,abs(TFdP_f))#*heaviside(-omega/2/pi+fisco,0.5)*heaviside(omega/2/pi-f_min,0.5))
# plt.loglog(freqdP_t,abs(TFdP_t))
# plt.show()

from matplotlib import rcParams,use,ticker
from numpy import pi
rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 600
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
#rcParams['figure.figsize'] = 15.5, 6.5
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 20
rcParams['font.size'] = 12
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 4
rcParams['legend.fontsize'] = 12
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] = 14

alpha,_,_=cav_parameters(cavity,K,l,Rmax,Rmin)
fig=plt.figure()
ax=plt.gca()

ax.set_box_aspect(1)
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax.plot(t_f,dP_f)
ax.plot(t,dP_t,linewidth=1,alpha=0.7)
ax.set_ylim(-1.2*max(abs(dP_t)),1.2*max(abs(dP_t)))
ax.set_title( cavity+' Cavity',fontdict={'fontweight':'bold'})
ax.set_xlabel("Time(s)")
ax.set_ylabel("Induced power (W)")

fig=plt.figure()
ax=plt.gca()
ax.set_box_aspect(1)
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax.loglog(omega/2/pi,abs(TFdP_f))
ax.loglog(freqdP_t,abs(TFdP_t),linewidth=1,alpha=0.7)
#ybottom = 10**-4;yup=10**3; ax.set_ylim(ybottom,yup)

# opt = dict(color='tab:blue', arrowstyle = 'simple,head_width=1.5,head_length=1.5',connectionstyle = 'arc3,rad=0',alpha=0.5)
# ax.annotate('',xy=(alpha[0]*c/l/2/pi, yup),xycoords='data',xytext=(alpha[0]*c/l/2/pi,ybottom),textcoords = 'data',arrowprops=opt,size=0.1)
# ax.annotate('',xy=(alpha[1]*c/l/2/pi, yup),xycoords='data',xytext=(alpha[1]*c/l/2/pi,ybottom),textcoords = 'data',arrowprops=opt,size=0.1)
# ax.annotate('',xy=(alpha[2]*c/l/2/pi, yup),xycoords='data',xytext=(alpha[2]*c/l/2/pi,ybottom),textcoords = 'data',arrowprops=opt,size=0.1)
# ax.annotate('',xy=(alpha[3]*c/l/2/pi, yup),xycoords='data',xytext=(alpha[3]*c/l/2/pi,ybottom),textcoords = 'data',arrowprops=opt,size=0.1)
# ax.annotate('',xy=(alpha[4]*c/l/2/pi, yup),xycoords='data',xytext=(alpha[4]*c/l/2/pi,ybottom),textcoords = 'data',arrowprops=opt,size=0.1)
#ax.set_xlim(10**6,3*10**8)
ax.set_title(cavity+' Cavity',fontdict={'fontweight':'bold'})
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel('Power spectrum (arb. unit)')

print(rmsP_t,rmsP_f)
plt.show()


