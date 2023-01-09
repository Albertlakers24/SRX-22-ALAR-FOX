import numpy as np
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto

#all the equations from Raymer

#Wing
#Surface, ft^2
S = m_mto * g / W_S_design
S_w = S/((ft_m)^2)
print(S_w)
#weight of fuel in wing, lb = 0
W_fw = 1
#aspect ratio
A = 12
#Sweep = 0
Lambda = 0
#dynamic pressure at cruise lb/ft^2
rho_cruise = ISA_calculator(280*FL_ft*ft_m,dt_cruise)[2]

#taper ratio

#t/c

#Ulitimate load = 1.5*limit load factor

#Design gross weight, lb (MTOW)
W_dg = m_mto/lbs_kg


W_wing = 0.036*S_w^(0.758)*W_fw^(0.0035)*(A/(np.cos(Lambda))^(2)^(0.6)*q^(0.006)*lambda^(0.04)*((100*tc)/np.cos(Lambda))^(-0.3)*(N_z*W_dg)^(0.49)

#Horizontal Tail

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


