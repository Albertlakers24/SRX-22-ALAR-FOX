import numpy as np
from Constants import *
from Initial_Aircraft_Sizing.Wing_planform import Sw, taper, c_mac, t_c_ratio
import matplotlib.pyplot as plt
from Airfoildata import Alpha_AF_xflr, CL_AF_xflr, Alpha_Wing, CL_Wing, Alpha_AF_xfoil, CL_AF_xfoil


eta_airfoil = 0.95          # Airfoil efficiency -- ADSEE Slides
CLmax_Clmax = 0.9          # Wing/Airfoil CL -- From graph - ADSEE II, Lecture 2, p.19
delta_alpha_CLMax = 25.46 * t_c_ratio     # for curved part of stall section -- From graph - p.21
mach = 0.45
Alpha0 = -4.45
Cl_max_AF = 1.61

def Calculate_beta(mach):
    beta = np.sqrt(1 - (mach ** 2))
    return beta

def Calculate_wingsweep(sweep_25, position, taper):
    sweep = np.degrees(np.arctan(np.tan(np.radians(sweep_25)) - (4/A*(position - 0.25) *(1-taper)/(1+taper))))
    return sweep

def Calculate_CL_alpha(beta, sweep):
    CL_alpha = (2*np.pi*A)/ (2 + np.sqrt(4 + ((A *beta/eta_airfoil)**2 * (1 + ((np.tan(np.radians(sweep)))**2)/(beta**2)))))
    return CL_alpha

def Calculate_CL_max(Cl_max, delta_Clmax):
    CL_max = CLmax_Clmax*Cl_max + delta_Clmax
    return CL_max

def Calculate_alpha_stall(CL_max, CL_alpha, alpha_0L):
    alpha_stall = (CL_max/CL_alpha) + alpha_0L + delta_alpha_CLMax
    return alpha_stall

def Calculate_Reynolds_number(rho, V, l, mu):
    Re = rho*V*l/mu
    return Re
mach_TO = 0.2
CL_max = Calculate_CL_max(Cl_max_AF,0)
beta_CR = Calculate_beta(mach)
beta_TO = Calculate_beta(mach_TO)
sweep_half = Calculate_wingsweep(0, 0.5, taper)
CL_alpha = Calculate_CL_alpha(beta_CR,sweep_half)
Alpha_stall = Calculate_alpha_stall(CL_max, CL_alpha,Alpha0 )
print('CL max', CL_max, 'CL_alpha', CL_alpha)
Alpha = np.arange(-10,23,0.25)
CL = np.array([])
for i in range(len(Alpha)):
    CL_value = CL_alpha * (np.radians(Alpha[i]) - np.radians(Alpha0))
    CL = np.append(CL, CL_value)

print('stall', np.degrees(Alpha_stall))
print('CL_max' , CL_max)

plt.plot(Alpha, CL, color = 'r', label = 'Wing DATCOM')
#plt.plot(Alpha_AF_xflr, CL_AF_xflr, color = 'g', label = 'Airfoil XFLR')
plt.plot(Alpha_AF_xfoil, CL_AF_xfoil, color = 'y', label = 'Airfoil XFOIL')
plt.plot(Alpha_Wing, CL_Wing, color = 'b', label = 'Wing XFLR5')
plt.xlabel('Alpha')
plt.ylabel('CL[-]')
plt.legend()
plt.show()

