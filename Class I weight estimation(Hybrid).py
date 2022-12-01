import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
R_norm = 1000 * 1852            #Range in meters
E = 45 * 60                     #Loiter endurance in seconds
V_cruise = 275 * 0.51444444     #m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
R_div = 0                       #m (TBD)  --> to be determined in literature
f_con = 5/100                   #-
e_kero = 42.9                   #MJ/kg Specific Energy Kerosene
e_atj = 43.2                    #MJ/kg Specific Energy SAF(ATJ)
e_lh2 = 142                     #MJ/kg Specific Energy Liquid Hydrogen
e_bat = 1.1                     #MJ/kg Specific Energy Battery (assuming 300Wh/kg)
s_takeoff_1524 = 1370           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 1370           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
density_0 =                     #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
density_1524=                   #1524m ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION

#Cdo calculations
Psi = 0.0075 #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97   #span efficiency factor (value based on Roelof reader p.46)
e = 1/((np.pi())*A*Psi+(1/phi))
Cfe =                           #Equivalent skin friction coefficient - depending on aircraft from empirical estimation
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe * Swet_S
CL = np.sqrt(np.pi()*Cd0*A*e)
CD = 2 * Cd0
#Aerodynamic Estimations
CL_to =
CD_to =
CL_0 =
CL_max =
# Stall Speed Requirement (CHECK WHICH RHO TO TAKE)
V_stall_land_1524 = np.sqrt(s_landing_1524 / 0.5915)
V_stall_takeoff_1524 = np.sqrt(s_takeoff_1524 / 0.5915)
W_S_takeoff = 1/2 * density_1524 * V_stall_takeoff_1524 **2 * CL_max_to
W_S_landing = 1/2 * density_1524 * V_stall_land_1524 **2 * CL_max_land

# Take off Distance Constraint
Rho = density_1524 / desnity_0 #
#W_P = TOP/ (W_S) * CL_to * Rho

# Rate of Climb Constraint

# Climb Gradent Constraint

# Cruise Speed Constraint

# Landing Distance Constraint
