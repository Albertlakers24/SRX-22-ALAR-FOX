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
rho_1524_rho0 = rho_1524/rho_0
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
eff_prop = 0.85              #Change with Literature
W_S = np.arange(1,3500,1)
##Cdo calculations
Psi = 0.0075                    #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #span efficiency factor (value based on Roelof reader p.46)
A = 12                           #Aspect Ratio (12-14) #Reference to ATR 72
e = 1/(np.pi*A*Psi+(1/phi))
Cfe = 0.0045                     #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1                     #(6.0-6.2) wetted area ratios -> depending on airframe structure
Cd0 = Cfe * Swet_S
#Aerodynamic Estimations
CL_max = 1.7                       #(1.2-1.8 twin engine) & (1.5-1.9 turboprop) max lift coefficient
CL_to = 1.9                        #Change with Estimate (1.7-2.1)
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
CL_land =  2.6                  #Change with Estimate (1.9-3.3)
TOP = 430                       #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC = 10.2                      #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0024                  #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_stall = 52                      #Change with CS25 or Requirement

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to * rho_1524_rho0 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_1524 * s_landing_1524/0.5847)/(2*0.95)
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_cruise_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_1524))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + CD_to/CL_to)*(np.sqrt((2/rho_1524)*(1/CL_to))))
#Stall Constraint
W_S_stall = 1/2 * rho_1524 * V_stall**2 * CL_max

plt.vlines(W_S_stall,0,0.5,'b',label="V_stall")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,0.5,'k',label ="Landing")
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,3500)
plt.ylim(0,0.4)
plt.xticks(np.arange(0,3501,500))
plt.yticks(np.arange(0,0.51,0.05))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.legend(loc = "upper right")
plt.grid()
plt.show()

'''#Mass Preliminary Calculation
W_P_design =
W_S_design = 
MTOW_design = 169800                    #N
P_max = MTOW_design / W_P_design
S = MTOW_design / W_S_design
#Degree of Hybridization of Energy (He) *Could be defined by each split point or total journey
t_cruise =
t_climb =
t_descent =
delta_v = V_cruise
E_total_climb = (MTOW_design*V_cruise)/ (Cl_cruise/Cd_cruise) * t_climb + (MTOW_design/g * delta_v**2) + MTOW_design*ROC*t_climb
E_total_cruise = (MTOW_design*V_cruise)/ (Cl_cruise/Cd_cruise) * t_cruise
E_total_descent = (MTOW_design*V_cruise)/ (Cl_cruise/Cd_cruise) * t_descent + (MTOW_design/g * delta_v**2) + MTOW_design*ROC*t_descent
E_total = E_total_climb+E_total_descent+E_total_cruise
E_nc =
He = E_nc / E_total         #Energy of non consumable(battery) / Total Energy
eta_stt = 0.85              #Efficiency chain from shaft-to-thrust
eta_btt = 0.95              #Efficiency chain from battery-to-thrust
NoD_ice = 2                 #Number of turboprop engines
NoD_em = 4                   #Number of electric motor engines
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
ddp = 0.8                   #Deep discharge protection
E_bat = 2.7*10**6           #Total Battery Energy per piece
m_fuel_ice = (1+tf)*P_ice*NoD_ice*BSFC*t_toal
m_bat = (1+ddp) * (E_nc/(eta_btt*E_bat))
m_OE = m_fuel_ice + m_bat + m_payload + 0.0009*MTOW_design**2 - 11.862*MTOW_design +49013           #Maximum Takeoff Mass'''
