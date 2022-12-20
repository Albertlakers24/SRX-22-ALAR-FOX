import numpy as np
#Constants

#1-NOx long CH4
#2-NOx long O3L
#3-CO2
Type = 2

#NOx
tau_CH = 12         #years
t = 13              #years


if Type==1:
    A = -5.16*10**-13
    E = 0.119*t
    RF = (A * np.exp(-t / tau_CH)) * (np.exp(t / tau_CH) * E - (tau_CH ** 2) * (np.exp(t / tau_CH) - np.exp(0)))
    print("RF_NOXlong_CH4=", RF)
if Type==2:
    A = -1.21*10**-13
    E = 0.119*t
    RF = (A * np.exp(-t / tau_CH)) * (np.exp(t / tau_CH) * E - (tau_CH ** 2) * (np.exp(t / tau_CH) - np.exp(0)))
    print("RF_NOXlong_O3L=", RF)
if Type==3:
