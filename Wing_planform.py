import numpy as np
import scipy as sp
from Class_I_Weight_Estimation import A, Swet_S, MTOW,V_cruise
from Lift_Drag_Polar import p, T, specific_gas_constant, b
#Switch for simple/double tapered wing
switch = 1                             # put 1 for simple tapered, 2 for double tapered

#Constants
gamma = 1.4                            # Specific heat ratio of gas
M_cross = 0.935                        # Technology factor for supercritical airfoils

#Calculation
Sw = 61                                # main wing area [m^2], change value base this on class I est.
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
    c_mac = (2/3)*c_r*((1+taper+taper**2)/1+taper)  # length of MAC
    y_mac = 0.5*(1/3)*(1+2*taper)/(1+taper)*b       # Spanwise location of MAC
    #Printing results
    print("t_c_ratio: ", t_c_ratio)
    print("Pressure [Pa]:", p)
    print("Cruise Mach number: ", M_cruise)
    print("Speed of sound at cruise", a_cruise)
    print("C_L_hat", C_Lhat)
    print("Cruise Speed [m/s]: ", V_cruise)
    print("cr",c_r)
    print("ct", c_t)

# add double taper wing later







