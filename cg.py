import numpy as np
import scipy as sp
from Fuselage import l_f
from Wing_planform import M_cruise, a_cruise, V_cruise

# x_cg = Sum(mass * distance)/Sum(distance)
#Switches
engine_placement = 1

# Some parameters
M_Dive = M_cruise + 0.09     # Diving mach number
VDive = a_cruise * M_Dive



#if statements
if engine_placement == 1:   # horizontal tailplane lever arm with the engine mounted at the wing
    l_H = 0.55 * l_f

else:
    l_H = 0.5 * l_f         # horizontal tailplane lever arm with the engine mounted at the wing


#Assumptions
#the fuselage is a cylindrical structure with uniform density

#masses
m_fuselage = 0.23 * np.sqrt(V_D*l_H/(w_f+h_f))*(S_f_wet**1.2)

#distances
x_fuselage = 0.5 * l_f


#
xpos=[0,]
ypos=[0,]

print(m_fuselage)