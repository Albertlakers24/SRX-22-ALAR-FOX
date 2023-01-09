import numpy as np
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *

#all the equations from Raymer

#Wing
#Surface, ft^2
S = m_mto * g / W_S_design
S_w = S/((ft_m)**2)

#weight of fuel in wing, lb = 0
W_fw = 1.
#aspect ratio
A = 12.
#Sweep = 0
Lambda = 0
#dynamic pressure at cruise lb/ft^2
rho_cruise = ISA_calculator(280*FL_ft*ft_m,dt_cruise)[2]
q_kgm2 = 0.5*rho_cruise*V_cruise**(2)
q_lbft2 = q_kgm2*23.73 #kg/m2 to lbs/ft2

#taper ratio
taper = 0.45

#t/c (assumption for now)
tc = 0.12

#Ulitimate load = 1.5*limit load factor (assumption for now)
Lim_load = 4.
N_z = 1.5*Lim_load

#Design gross weight, lb (MTOW)
W_dg = m_mto/lbs_kg


W_wing_lb = 0.036*S_w**(0.758)*W_fw**(0.0035)*(A/(np.cos(Lambda))**(2))**(0.6)*q_lbft2**(0.006)*taper**(0.04)*((100*tc)/np.cos(Lambda))**(-0.3)*(N_z*W_dg)**(0.49)
W_wing = W_wing_lb * lbs_kg
print('Weight Wing kg =', W_wing)


#Horizontal Tail
from Initial_Aircraft_Sizing.
#S_ht
#tc_ht
#taper_ht
taper_ht =
W_hor_tail = 0.016*(N_z*Wdg)**(0.414)*q_lbft2**(0.168)*S_ht*(0.896)*((100*tc_ht)/np.cos(Lambda))**(-0.12)*(A/(np.cos(Lambda))**(2))**(0.043)*taper_ht**(-0.02)
print('Weight Horizontal Tail = ', W_hor_tail)

#Vertical Tail

#Fuselage

#Main Landing Gear

#Nose Landing Gear

#Installed Engine
#(From our own estimates and NOT Raymer)

#Fuel System
#(From our own estimates and NOT Raymer)

#Flight Controls

#Hydrolics

#Electrical

#Avionics

#Airconditioning and Anti-Icing

#Furnishing


