from cProfile import label
from matplotlib import rcParams,ticker
from matplotlib.pyplot import *
from numpy import fft,pi
import pickle

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
rcParams['text.usetex'] = True


with open('./SIMTIMEM5D25TEM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM5D25TEM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(1)
ax1=subplot(1,2,1)
plot(t,dP_t,'b')
axis([0,6e-5,-5e-9,5e-9])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
axis([1e4,1e9,1e-42,1e-15])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
suptitle(r'\textbf{TEM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh5TEM.png")

with open('./SIMTIMEM5D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM5D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(2)
ax1=subplot(1,2,1)
plot(t,dP_t,'b')
axis([0,6e-5,-5e-9,5e-9])
xlabel(r'Time(s)')
ylabel(r'Induced EM Power (W)')
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax1.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True, useOffset=False))
ax2=subplot(1,2,2)
loglog(omega/2/pi,abs(TFdP_f)**2,'g')
axis([1e4,1e9,1e-42,1e-15])
xlabel(r'Frequency (Hz)')
ylabel(r' $\left\vert \mathcal{P} (\omega) \right\vert^2 (W^2)$')
ax2.set_box_aspect(1)
suptitle(r'\textbf{TM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh5TM.png")


fig=figure(3)
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
suptitle(r'\textbf{TM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
handles,labels=ax1.get_legend_handles_labels()
fig.legend(handles,labels,loc='upper right')
savefig("pbh5comp.png")

with open('./SIMTIMEM3D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM3D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(4)
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
suptitle(r'\textbf{TM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh3TM.png")


with open('./SIMTIMEM7D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,t,dP_t,rmsP_t = pickle.load(f)
    f.close()
with open('./SIMFREQM7D25TM.pkl','rb') as f:  # Python 3: open(..., 'rb')
    M,omega,TFdP_f,rmsP_f,t_f,dP_f = pickle.load(f)
    f.close()

figure(5)
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
suptitle(r'\textbf{TM Cavity}')
subplots_adjust(left=0.07,bottom=0.11,right=0.97,top=0.933,wspace=0.236,hspace=0.2)
savefig("pbh7TM.png")




