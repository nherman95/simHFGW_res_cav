from math import log2,pi
from cavities import *
from scipy.integrate import quadrature,romberg,simps,romb,trapz
from scipy.stats import norm
import matplotlib.pyplot as plt
npts=1e7
t=linspace(0,1000,2**round(log2(npts))+1)
print(sqrt(romb(cos(t)**2,t[1]-t[0])/1000))
print(1/sqrt(2))
print(sqrt(sum(cos(t)**2)/(2**round(log2(npts))+1)))
tf=fft.rfft(cos(t))
freq=fft.rfftfreq(len(t))
print(sum((abs(tf)/len(t))**2))
print(freq)