import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *
from Initial_Aircraft_Sizing.Wing_planform import b

# Density
rho = 0
density = 0
if density == 0:
    rho = rho_0 # density at sea level
elif density == 1:
    rho = rho_cruise
else:
    rho = rho_loiter # add loiter density to Constants file

# MANEUVER DIAGRAM DESIGN
#Max lift coefficient
if rho == 0:
    Clmax = CL_max_landing
elif rho == 1:
    Clmax = CL_max_cruise
else:
    Clmax = CL_max_loiter # add loiter Clmax to Constants file

#Load factor
n = 2.1 + (24000 / (m_mto + 10000))

if n < 2.5:
    n_max = 2.5
if n > 3.8:
    n_max = 3.8

n_min = -1

#Cruise speed EAS
V_C = V_cruise
V_C = V_C * np.sqrt(rho/rho_0)

#Dive speed EAS
V_D1 = V_C/0.8
V_D2 = #Transform V_C to mach number, add 0.05, and transform back to m/s
V_D = min(V_D1, V_D2)
V_D = V_D * np.sqrt(rho/rho_0)

#Stall speed
V_s = np.sqrt((2 * 1 * W_S_design) / (rho_0 * Clmax)) * np.sqrt(rho / rho_0)

