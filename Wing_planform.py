import numpy as np
#import scipy as sp
from Lift_Drag_Polar import p, T, specific_gas_constant
#from Fuselage import prop_choice
#Switch for simple/double tapered wing
switch = 1                            # put 1 for simple tapered, 2 for double tapered
prop_type = 4                        # 1 = H2 combustion, 2 = fuel cell, 3 = EH Series, 4 = EH Parallel Sereis
#Constants
gamma = 1.4                            # Specific heat ratio of gas
M_cross = 0.935                        # Technology factor for supercritical airfoils
g = 9.80665
# Parameters per aircraft
if prop_type == 1:
    MTOM = 19200    # LH2 combustion in kg
if prop_type == 2:
    MTOM = 19400    # hydrogen fuel cell
if prop_type == 3:
    MTOM = 25300    # battery fuel hybrid in series
if prop_type == 4:
    MTOM = 24100    # battery fuel hybrid in parallel series

Sw = MTOM *g / 3171     # To be updated!!
A = 12
b = np.sqrt(A*Sw)
print('Sw', Sw, 'b', b)

# Mission charecteristics
V_cruise = 141.471   # in m/s
#Calculation
#Sw = 61                               # main wing area [m^2], change value base this on class I est.
a_cruise = np.sqrt(gamma*specific_gas_constant*T)          # Speed of sound at cruise altitude  [m/s]
M_cruise = V_cruise / a_cruise         # Mach number during cruise
M_dd = M_cruise + 0.03                 # Drag-divergence Mach number
taper = 0.45                           # Taper ratio (from 0 to 1), optimum taper is 0.45 for unswept
                                       # wing to achieve elliptical lift dist. check Raymer

# for simple tapered wing, used for deciding on the airfoil
if switch == 1:
    c_r = 2*Sw/((1+taper)*b)           # Tip chord  [m]
    c_t = c_r * taper                  # Root chord [m]
    qhat = 0.5 * gamma * p * (M_cruise**2)
    C_Lhat = MTOM/(qhat*Sw)            # Cruise lift coefficient
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
    print("Spanwise position of MAC", y_mac)

if switch == 2:                         # For double tapered wing

    # c_r = 2*Sw/((1+taper)*b)
    eta_k = 0.4                         # relative span position of kink, need to determine this value later
                                        # Depending on the postition of the propellor
    y_k = b * eta_k / 2                 # Spanwise position of the kink
    taper_inner = 1                     # Taper ratio of the inner wing [to be changed]
    taper_outer = 0.5                   # Taper ratio of the outer wing [to be changed]
    c_r = (2*Sw/b) / ((eta_k*(1-taper_outer))+taper_inner+taper_outer)
    c_k = c_r * taper_inner             # chord length at kink
    c_t = c_k * taper_outer             # chord length at tip
    c_mac_inner = 2/3*c_r*(1+taper_inner+taper_inner**2)/(1+taper_inner)        # MAC of inner wing [m]
    c_mac_outer = 2/3*c_r*(1+taper_outer+taper_outer**2)/(1+taper_outer)        # MAC of outer wing [m]
    Sw_inner = b/2*eta_k*(c_r)*(1+taper_inner)                                  # inner wing area [m]
    Sw_outer = Sw - Sw_inner                                                    # outer wing area [m]
    c_mac = (c_mac_inner * Sw_inner + c_mac_outer * Sw_outer) / Sw              # MAC of wing [m]
    y_mac_inner = (b/2*eta_k)/3*((1+2*taper_inner)/(1+taper_inner))             # spanwise location of MAC of inner wing [m]
    y_mac_outer = (b/2*(1-eta_k))/3*((1+2*taper_outer)/(1+taper_outer))+(b/2*eta_k) # spanwise location of MAC of outer wing [m]
    y_mac = (y_mac_inner * Sw_inner + y_mac_outer * Sw_outer) / Sw              # spanwise location of MAC of wing [m]


    # Printing results
    print("root chord length [m]", c_r)
    print("kink chord length [m]", c_k)
    print("tip chord length  [m]", c_t)
    print("MAC of outer wing [m]", c_mac_outer)
    print("MAC of inner wing [m]", c_mac_inner)
    print("Inner wing area  [m]", Sw_inner)
    print("Outer wing area  [m]", Sw_outer)
    print("MAC length  [m]", c_mac)
    print("spanwise position of inner wing MAC  [m]", y_mac_inner)
    print("spanwise position of outer wing MAC  [m]", y_mac_outer)
    print("spanwise position of total wing MAC  [m]", y_mac)









