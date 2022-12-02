import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
R_norm = 1000 * 1852            #Range in meters
E = 45 * 60                     #Loiter endurance in seconds
V_cruise = 275 * 0.51444444     #knt -> m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
R_div = 0                       #m (TBD)  --> to be determined in literature
f_con = 5/100                   #-
e_kero = 42.9                   #MJ/kg Specific Energy Kerosene
e_atj = 43.2                    #MJ/kg Specific Energy SAF(ATJ)
e_lh2 = 142                     #MJ/kg Specific Energy Liquid Hydrogen
e_bat = 1.1                     #MJ/kg Specific Energy Battery (assuming 300Wh/kg)
s_takeoff_1524 = 1370           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 1370           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
rho_0=                          #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524=                       #1524m ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
W_S = np.arange(1,3000,1)
##Cdo calculations
Psi = 0.0075 #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97   #span efficiency factor (value based on Roelof reader p.46)
e = 1/((np.pi())*A*Psi+(1/phi))
Cfe =                           #Equivalent skin friction coefficient
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe * Swet_S
CL = np.sqrt(np.pi()*Cd0*A*e)
CD = 2 * Cd0
#Aerodynamic Estimations
CL_to =
CD_to =
CL_0 =
CL_max =
## THE DOCUMENT HAS P/W INSTEAD OF W/P
# Stall Speed Requirement (CHECK WHICH RHO TO TAKE)
V_stall_land_1524 = np.sqrt(s_landing_1524 / 0.5915)
V_stall_takeoff_1524 = np.sqrt(s_takeoff_1524 / 0.5915)
stall_criteria_takeoff = 1/2 * density_1524 * V_stall_takeoff_1524 **2 * CL_max_to
stall_criteria_landing = 1/2 * density_1524 * V_stall_land_1524 **2 * CL_max_land

# Take off Distance Constraint
Rho = density_1524 / density_0 #
W_P_TOP = TOP/ ((W_S) * CL_to * Rho)

# Rate of Climb Constraint

# Climb Gradent Constraint
W_P_cv = eff_prop / (np.sqrt(W_S)*(climb_grad + CD_to/CL_to)*(np.sqrt((2/density_1524)*(1/CL_to))))
# Cruise Speed Constraint
W_P_cru
# Landing Distance Constraint
W_P_land
#Degree of Hybridization of Power(Hp)
# Choice between Parallel and Series needs to be made
#If Parallel:
H_p_para = P_em_max / P_max
#If Series:
H_p_ser = P_em_max / P_ice_max

#Degree of Hybridization of Energy (He) *Could be defined by each split point or total journey
He = E_nc / E_total         #Energy of non consumable(battery) / Total Energy
#Create Wing Loading diagram and split points
plt.vlines(stall_criteria_takeoff,0,1,'b',label = "V_stall Takeoff")
plt.plot(W_S,1/W_P_TOP,'r',label = "Takeoff Constraint")
plt.plot(W_S,1/W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,1/W_P_cv,'g',label = "Climb Gradient Constraint")
plt.plot(W_S,1/W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,1/W_P_land,'y',label = "Landing Constraint")
plt.xlabel("W/S")
plt.ylabel("W/P")
plt.ylim(0,1)
plt.legend(loc="upper left")
plt.show()
#Create split point by ratio
