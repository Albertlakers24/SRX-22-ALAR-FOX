import numpy as np
from Constants import *
from Initial_Aircraft_Sizing.Wing_planform import Sw, taper, c_mac
import matplotlib.pyplot as plt

eta_airfoil = 0.95          # Airfoil efficiency -- ADSEE Slides
CLmax_Clmax = 0.88          # Wing/Airfoil CL -- From graph - ADSEE II, Lecture 2, p.19
delta_alpha_CLMax = 2.1     # for curved part of stall section -- From graph - p.21
delta_Clmax_TO = 0.01       # compressibility -- from graph - based on mach number - p.20
delta_Clmax_land = 0.04     # compressibility -- from graph - based on mach number - p.20

def Calculate_beta(mach):
    beta = np.sqrt(1 - (mach ** 2))
    return beta

def Calculate_wingsweep(sweep_25, position, taper):
    sweep = np.degrees(np.arctan(np.tan(np.radians(sweep_25)) - (4/A*(position - 0.25) *(1-taper)/(1+taper))))
    return sweep

def Calculate_CL_alpha(beta, sweep):
    CL_alpha = 2*np.pi*A/ (2 +np.sqrt((A *beta/eta_airfoil)**2 * (1 + (np.tan(np.radians(sweep)))**2/(beta**2))))
    return CL_alpha

def Calculate_CL_max(Cl_max, delta_Clmax):
    CL_max = CLmax_Clmax*Cl_max + delta_Clmax
    return CL_max

def Calculate_alpha_stall(CL_max, CL_alpha, alpha_0L):
    alpha_stall = CL_max/CL_alpha + alpha_0L + delta_alpha_CLMax
    return alpha_stall

def Calculate_Reynolds_number(rho, V, l, mu):
    Re = rho*V*l/mu
    return Re

sweep_LE = Calculate_wingsweep(0, 0)
print(sweep_LE)
CL_max = Calculate_CL_max(1.61,-0.45)
print(CL_max)
mach = 0.45
beta = Calculate_beta(mach)
sweep_half = Calculate_wingsweep(0, 0.5)
CL_alpha = Calculate_CL_alpha(beta,sweep_half)
CL = []
Alpha = np.arange(-10,30,0.25)
Alpha0 = -4.5
for i in range(len(Alpha)):
    CL_value = CL_alpha * (np.radians(Alpha[i]) - np.radians(Alpha0))
    CL.append(CL_value)

print(CL, Alpha)
print(np.where(CL == 0), Alpha[np.where(CL==0)])
plt.plot(Alpha, CL)
plt.show()
print(sweep_half)