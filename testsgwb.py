from cavities import Sh_toy,freq_appoach
from numpy import pi,linspace,sqrt
from matplotlib import rcParams
rcParams["font.family"] = "serif"
rcParams["font.weight"] = 500
rcParams["font.serif"] = 'STIXgeneral'
rcParams["mathtext.fontset"] = "stix"
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.1

rcParams['figure.autolayout'] = False
rcParams['figure.figsize'] = 9.5, 8
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 36
rcParams['font.size'] = 16
rcParams['lines.linewidth'] = 2.5
rcParams['lines.markersize'] = 0.5
rcParams['legend.fontsize'] = 14
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22

l=1
cavity='TEM'
Rmax=5
K=5
omegaGW=1e-10
fmax =1e10
fmin = 1e5
dt = 0.25*fmax**-1
nbrpt = 1/fmin/dt
print(nbrpt)
freq = linspace(fmin,fmax,int(nbrpt))
omega = 2*pi*freq
fcut=1e8
test=sqrt(2)*Sh_toy(freq,omegaGW,fcut)
msP_f,TFdP_f_p,dP_f,tdE_f=freq_appoach(test,omega,cavity,Rmax)