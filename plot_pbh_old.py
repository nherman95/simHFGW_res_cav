from cProfile import label
from matplotlib import rcParams,ticker
from matplotlib.pyplot import *
from numpy import fft,pi
from math import log
import pickle
from cavities import RMS_compute,fchirp, gw_binary_time_signal,gw_binary_freq_signal

rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 400
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 16, 9
rcParams['axes.labelsize'] = 32
rcParams['figure.titlesize'] = 42
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 8
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
#rcParams['text.usetex'] = True


with open('./SIMTIMEM5D25TEM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM5D25TEM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()
fsize=28
figure(1)
for i, label in enumerate(('a)', 'b)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
print(rmsP_t,rmsP_f,RMS_compute(fchirp,2200/25/M,2200/M,'TEM',5,M))
ax1=subplot(1,2,1)
plot(t,dP_t,'b')
axis([0,6e-5,-5e-9,5e-9])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
annotate(r"resonance triggered",(5.1e-5,5e-10), xytext=(1e-5,2e-9),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=fsize,fontweight='bold')
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
axis([1e4,1e9,1e-42,1e-15])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
suptitle(r'TEM Cavity',fontweight='bold')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh5TEM.pdf")

with open('./SIMTIMEM5D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM5D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(2)
for i, label in enumerate(('c)', 'd)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
ax1=subplot(1,2,1)
plot(t,dP_t,'b')
axis([0,6e-5,-5e-9,5e-9])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
annotate(r"resonance triggered",(5.1e-5,5e-10), xytext=(1e-5,2e-9),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=fsize, fontweight='bold')
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
axis([1e4,1e9,1e-42,1e-15])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
suptitle(r'TM Cavity', fontweight='bold')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh5TM.pdf")


fig=figure(3)
for i, label in enumerate(('a)', 'b)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
ax1=subplot(1,2,1)
TFdP_t = fft.rfft(dP_t,n=len(dP_t))*2/len(dP_t)
freqdP_t = fft.rfftfreq(len(dP_t),t[1]-t[0])
plot(t,dP_t,'b',label='Time domain code')
plot(t_f,dP_f,'g',alpha=0.6,label='Frequency domain code')
axis([0,6e-5,-5e-9,5e-9])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
loglog(freqdP_t,abs(TFdP_t)**2,'b',alpha=0.6)
axis([1e4,1e9,1e-42,1e-15])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
suptitle(r'TM Cavity', fontweight='bold')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
handles,labels=ax1.get_legend_handles_labels()
fig.legend(handles,labels,loc='upper right')
savefig("pbh5comp.pdf")

with open('./SIMTIMEM3D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM3D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(4)
for i, label in enumerate(('a)', 'b)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.2, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
print(rmsP_t,rmsP_f,RMS_compute(fchirp,2200/25/M,2200/M,'TM',5,M))
ax1=subplot(1,2,1)
plot(t,dP_t,'b')
axis([0,6e-3,-6e-12,6e-12])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
axis([1e1,1e9,1e-52,1e-28])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
annotate(r"$\propto \omega^{11/6}$",(5e5,1e-32), xytext=(1e3,1e-31),arrowprops=dict(facecolor='black', shrink=0.05),fontsize=fsize,fontweight='bold')
suptitle(r'TM Cavity', fontweight='bold')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh3TM.pdf")


with open('./SIMTIMEM7D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM7D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(5)
for i, label in enumerate(('a)', 'b)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.15, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
print(rmsP_t,rmsP_f,RMS_compute(fchirp,2200/25/M,2200/M,'TEM',5,M))
ax1=subplot(1,2,1)
plot(t,dP_t,'b')
axis([0,6e-7,-2e-11,2e-11])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
axis([1e6,1e11,1e-52,1e-25])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
suptitle(r'TM Cavity', fontweight='bold')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh7TM.pdf")




fig=figure(6)
for i, label in enumerate(('a)', 'b)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.2, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
M=1e-5
fisco=2200/M
hp,hx,t=gw_binary_time_signal(M,1e9,fisco/25,1.0/2**(int(round(log(fisco,2)))+2))
ax1=subplot(1,2,1)
TFh_t = fft.rfft(hp,n=len(hp))*2/len(hp)
freqdP_t = fft.rfftfreq(len(hp),t[1]-t[0])
plot(t,hp,'b',label='Time domain code')
#plot(t_f,dP_f,'g',alpha=0.6,label='Frequency domain code')
axis([0,6e-5,-1.5e-28,1.5e-28])
xlabel(r'Time(s)')
ylabel(r'Strain')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
#loglog(omega/2/pi,abs(TFdP_f)**2,'g')
loglog(freqdP_t,abs(TFh_t)**2,'b')
axis([1e4,1e9,1e-67,1e-59])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \tilde{h} (\omega) \right\vert^2 $')
ax2.set_box_aspect(1)
#suptitle(r'\textbf{TM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.286,hspace=0.2)
savefig("ht.pdf")

fig=figure(7)
for i, label in enumerate(('c)', 'd)')):
    ax = gcf().add_subplot(1,2,i+1)
    ax.text(-0.2, 1.05, label, transform=ax.transAxes,
      fontsize=24, fontweight='bold', va='top')
M=1e-5
fisco=2200/M
TFhp,omega=gw_binary_freq_signal(M,1e9,fisco/25,fisco,1/t[-1])
ax1=subplot(1,2,1)
h = fft.irfft(TFhp,n=len(hp))*len(TFhp)
plot(t,h,'g')
#plot(t_f,dP_f,'g',alpha=0.6,label='Frequency domain code')
axis([0,6e-5,-1.5e-28,1.5e-28])
xlabel(r'Time(s)')
ylabel(r'Strain')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFhp)**2,'g')
#loglog(freqdP_t,abs(TFh_t)**2,'b')
axis([1e4,1e9,1e-67,1e-59])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \tilde{h} (\omega) \right\vert^2 $')
ax2.set_box_aspect(1)
#suptitle(r'\textbf{TM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.286,hspace=0.2)
savefig("hf.pdf")
# handles,labels=ax1.get_legend_handles_labels()
# fig.legend(handles,labels,loc='upper right')
