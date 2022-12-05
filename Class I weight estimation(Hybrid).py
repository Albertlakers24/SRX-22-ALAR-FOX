import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
lambda_trop = -6.5/1000
R = 287.0528
R_norm = 1000 * 1852            #Range in meters
E = 45 * 60                     #Loiter endurance in seconds
V_cruise = 275 * 0.51444444     #knt -> m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
R_div = 0                       #m (TBD)  --> to be determined in literature
f_con = 5/100                   #-
s_takeoff_1524 = 1370           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 1370           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
rho_0= 1.225                    #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524= 1.01893               #1524m ISA + 10 ◦C day (kg/m3)
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
eff_prop = 0.8              #Change with Literature
W_S = np.arange(1,3000,1)
##Cdo calculations
Psi = 0.0075                    #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #span efficiency factor (value based on Roelof reader p.46)
A =                             #Aspect Ratio
e = 1/(np.pi*A*Psi+(1/phi))
Cfe =                           #Equivalent skin friction coefficient
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe * Swet_S
CL = np.sqrt(np.pi()*Cd0*A*e)
#Aerodynamic Estimations
CL_max =                        #Change with Literature
CL_to =                         #Change with Estimate
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
CL_land =                       #Change with Estimate
TOP =                           #Change with Literature Reference to slide
ROC =                           #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V =                         #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_stall =                       #Change with CS25 or Requirement
W_S = np.arange(1,1500,1)
#Stall Constraint
W_S_stall = 1/2 * rho_0 * V_stall**2 * CL_max
#Takeoff Constraint
W_P_TOP = TOP/ (W_S) * CL_to #* rho_rho0 (Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_0 * s_landing_1524/0.5915)/(2*0.95) #Change to CS25 regulation
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_cruise_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + CD_to/CL_to)*(np.sqrt((2/rho_0)*(1/CL_to))))
plt.vlines(W_S_stall,0,0.4,'b',label="V_stall")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,0.4,'k',label ="Landing")
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,1500)
plt.ylim(0,0.4)
plt.xticks(np.arange(0,1501,500))
plt.yticks(np.arange(0,0.41,0.05))
plt.legend(loc = "upper right")
plt.grid()
plt.show()
#Degree of Hybridization of Power(Hp)
# Choice between Parallel and Series needs to be made
#If Parallel:
H_p_para = P_em_max / P_max
#If Series:
H_p_ser = P_em_max / P_ice_max

#Degree of Hybridization of Energy (He) *Could be defined by each split point or total journey
He = E_nc / E_total         #Energy of non consumable(battery) / Total Energy

#Mass Preliminary Calculation
P_max = MTOW_desing / W_P_design
S = MTOW / W_S_design
