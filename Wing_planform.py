import numpy as np
import scipy as sp
from Class_I_Weight_Estimation import b, A, Swet_S, MTOW,V_cruise
from Lift_Drag_Polar import p, T, specific_gas_constant

#Constants
gamma = 1.4                            # Specific heat ratio of gas
M_cross = 0.935                        # Technology factor for supercritical airfoils

#Calculation
a_cruise = np.sqrt(gamma*specific_gas_constant*T)          # Speed of sound at cruise altitude  [m/s]
M_cruise = V_cruise / a_cruise         # Mach number during cruise
M_dd = M_cruise + 0.03                 # Drag-divergence Mach number
taper = 0.4                            # Taper ratio (from 0 to 1), 0.4 for unswept wing. check Raymer
c_r = 2*Swet_S/((1+taper)*b)           # Tip chord  [m]
c_t = c_r * taper                      # Root chord [m]
qhat = 0.5 * gamma * p * (M_cruise**2)
C_Lhat = MTOW/(qhat*Swet_S)           # Cruise lift coefficient
t_c_ratio = min(0.18, ((M_cross-M_dd)-0.115*(C_Lhat**1.5))) # thickness to chord ratio


#printing results
print("t_c_ratio: ", t_c_ratio)
print("Pressure [Pa]:", p)
print("Cruise Mach number: ", M_cruise)
print("Speed of sound at cruise", a_cruise)
print("C_L_hat", qhat)
print("Cruise Speed [m/s]: ", V_cruise)
