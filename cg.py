import numpy as np
import scipy as sp
from Fuselage import l_f, D_outer, S_f_wet, Mass_fuselage
from Wing_planform import M_cruise, a_cruise, V_cruise

# x_cg = Sum(mass * distance)/Sum(distance)
#Switches
engine_placement = 1

# Some parameters
M_Dive = M_cruise + 0.09     # Diving mach number
VDive = a_cruise * M_Dive
w_f = D_outer                # Width of fuselage
h_f = w_f                    # Height of fuselage, assumed to be a circular fuselage



#if statements
if engine_placement == 1:   # horizontal tailplane lever arm with the engine mounted at the wing
    l_H = 0.55 * l_f
else:
    l_H = 0.5 * l_f         # horizontal tailplane lever arm with the engine mounted at the wing




#Assumptions
#the fuselage is a cylindrical structure with uniform density

#masses


#distances
x_fuselage = 0.5 * l_f


#
xlst=[0,]
mlst=[0,Mass_fuselage]

print(m_fuselage)
