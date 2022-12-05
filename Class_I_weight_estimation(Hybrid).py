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
R_div = 200*1000                #m (TBD)  --> to be determined in literature
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
A =                             #Aspect Ratio (12-14)
e = 1/(np.pi*A*Psi+(1/phi))
Cfe =                           #Equivalent skin friction coefficient
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe * Swet_S
#Aerodynamic Estimations
CL_max =                        #(1.2-1.8 twin engine) & (1.5-1.9 turboprop) max lift coefficient
CL_to =                         #Change with Estimate (1.9-3.3)
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
CL_land =                       #Change with Estimate (1.7-2.1)
TOP =                           #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC =                           #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V =                         #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_stall =                       #Change with CS25 or Requirement
W_S = np.arange(1,1500,1)

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to * rho_rho0 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_0 * s_landing/0.5915)/(2*0.95) #Change to CS25 regulation
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_cruise_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + CD_to/CL_to)*(np.sqrt((2/rho_0)*(1/CL_to))))
#Stall Constraint
W_S_stall = 1/2 * rho_0 * V_stall**2 * CL_max

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

#Mass Preliminary Calculation
MTOW_design=
W_P_design =
W_S_design =
P_max = MTOW_design / W_P_design
S = MTOW_design / W_S_design
#Degree of Hybridization of Energy (He) *Could be defined by each split point or total journey
t_toal =
t_cruise =
delta_v = 
E_nc =
E_total = (MTOW_design*V_cruise)/ (Cl_cruise/Cd_cruise) * t_cruise + (MTOW_design/g * delta_v**2) + MTOW_design*ROC*(t_total - t_cruise)
He = E_nc / E_total         #Energy of non consumable(battery) / Total Energy
eta_stt =                   #Efficiency chain from shaft-to-thrust
eta_btt =                   #Efficiency chain from battery-to-thrust
NoD_ice =                   #Number of turboprop engines
NoD_em =                    #Number of electric motor engines
P_ice = (E_total - E_nc)/ (eta_stt * t_toal * NoD_ice)
P_em = E_nc/ (eta_btt * t_toal * NoD_em)
#Degree of Hybridization of Power(Hp)
# Choice between Parallel and Series needs to be made
#If Parallel:
H_p_para = P_em_max / P_max
#If Series:
H_p_ser = P_em_max / P_ice_max
tf =                        #Trap fuel time step
BSFC=                       #Brake-specific fuel consumption
ddp =                       #Deep discharge protection
E_bat =                     #Total Battery Energy
m_fuel_ice = (1+tf)*P_ice*NoD_ice*BSFC*t_toal
m_bat = (1+ddp) * (E_nc/(eta_btt*E_bat))
m_OE = m_fuel_ice + m_bat + m_payload + 0.0009*MTOW_design**2 - 11.862*MTOW_design +49013           #Maximum Takeoff Mass
