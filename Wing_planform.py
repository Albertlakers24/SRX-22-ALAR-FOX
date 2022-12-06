import numpy as np
import scipy as sp
from Class_I_Weight_Estimation import MTOW,V_cruise
from Lift_Drag_Polar import p, T, specific_gas_constant, b
#Switch for simple/double tapered wing
switch = 2                             # put 1 for simple tapered, 2 for double tapered

#Constants
gamma = 1.4                            # Specific heat ratio of gas
M_cross = 0.935                        # Technology factor for supercritical airfoils

#Calculation
Sw = 61                               # main wing area [m^2], change value base this on class I est.
a_cruise = np.sqrt(gamma*specific_gas_constant*T)          # Speed of sound at cruise altitude  [m/s]
M_cruise = V_cruise / a_cruise         # Mach number during cruise
M_dd = M_cruise + 0.03                 # Drag-divergence Mach number
taper = 0.45                           # Taper ratio (from 0 to 1), 0.45 for unswept wing. check Raymer


# for simple tapered wing, used for deciding on the airfoil
if switch == 1:
    c_r = 2*Sw/((1+taper)*b)           # Tip chord  [m]
    c_t = c_r * taper                  # Root chord [m]
    qhat = 0.5 * gamma * p * (M_cruise**2)
    C_Lhat = MTOW/(qhat*Sw)            # Cruise lift coefficient
    t_c_ratio = min(0.18, ((M_cross-M_dd)-0.115*(C_Lhat**1.5))) # thickness to chord ratio
    c_mac = (2/3)*c_r*((1+taper+taper**2)/(1+taper))  # length of MAC
    y_mac = 0.5*(1/3)*(1+2*taper)/(1+taper)*b       # Spanwise location of MAC
    #Printing results
    print("t_c_ratio: ", t_c_ratio)
    print("Pressure [Pa]:", p)
    print("Cruise Mach number: ", M_cruise)
    print("Speed of sound at cruise", a_cruise)
    print("C_L_hat", C_Lhat)
    print("Cruise Speed [m/s]: ", V_cruise)
    print("MAC: ", c_mac)
    print("root chord: ",c_r)
    print("tip chord: ", c_t)

if switch == 2:                         # For double tapered wing

    # c_r = 2*Sw/((1+taper)*b)
    eta_k = 0.4                         # relative span position of kink, need to determine this value later
    y_k = b * eta_k / 2                 # Spanwise position of the kink
    taper_inner = 1                     # Taper ratio of the inner wing [to be changed]
    taper_outer = 0.5                   # Taper ratio of the outer wing [to be changed]
    c_r = (2*Sw/b) / ((eta_k*(1-taper_outer))+taper_inner+taper_outer)
    c_k = c_r * taper_inner             # chord length at kink
    c_t = c_k * taper_outer             # chord length at tip
    #add MAC location
    #add MAC length
    # Printing results
    print("root chord length [m]", c_r)
    print("kink chord length [m]", c_k)
    print("tip chord length  [m]", c_t)








