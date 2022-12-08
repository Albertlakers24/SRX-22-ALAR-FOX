import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
lambda_trop = -6.5/1000
R = 287.0528
V_cruise = 280 * 0.51444444     #knt -> m/s (TAS)
h_cruise = 200*100 * 0.3048     #m
s_takeoff_1524 = 1315           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 1315           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
rho_0= 1.225                    #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524= 1.01893               #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
eff_prop = 0.85                 #Change with Literature
m_payload = 8124
W_S = np.arange(1,5000,1)
##Cdo calculations
Psi = 0.0075                    #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #span efficiency factor (value based on Roelof reader p.46)
A = 12                           #Aspect Ratio (12-14) #Reference to ATR 72
e = 1/(np.pi*A*Psi+(1/phi))
Cfe = 0.0045                     #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1                     #(6.0-6.2) wetted area ratios -> depending on airframe structure
Cd0 = Cfe * Swet_S
#Aerodynamic Estimations
CL_max = 1.9                       #(1.2-1.8 twin engine) & (1.5-1.9 turboprop) max lift coefficient
CL_to = 2.1                        #Change with Estimate (1.7-2.1)
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
CL_land = 2.6                      #Change with Estimate (1.9-3.3)
CD_land = Cd0 + (CL_to**2 /(np.pi * A* e))
TOP = 500                          #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC = 6.88                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach = 141* 0.514444         #Change with CS25 or Requirement
a = 0.5088
b = 1199.7

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to #* rho_1524_rho0 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_0 * s_landing_1524/0.5847)/(2*0.95)
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_cruise_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + (CD_to/CL_to))*(np.sqrt((2/rho_0)*(1/CL_to))))
#Stall Constraint
W_S_approach = 1/2 * rho_0 * V_approach**2 * CL_land

plt.vlines(W_S_approach,0,100,'b',label="V_approach")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,100,'k',label ="Landing")
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,5000)
plt.ylim(0,0.5)
plt.xticks(np.arange(0,5001,500))
plt.yticks(np.arange(0,0.5,0.05))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.legend(loc = "upper right")
plt.grid()
plt.show()

#Mass Preliminary Calculation
W_P_design = 0.053
W_S_design = 3767
m_turboprop = 1074.5/2
m_em = 50
m_propeller = 17
MTOW_design = 274732                   #N
P_max = MTOW_design / W_P_design
S = MTOW_design / W_S_design
print(P_max/1000, "kW Max Power")
print(S,"m^2 Surface Area ")

#Energy Calculation
t_cruise_500 = 4476
t_cruise_full = 17568
t_climb1 = 222
t_climb2 = 600
t_climb3 = 975
t_descent1 = 503
t_descent2 = 377
t_total_500 = t_cruise_500 + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
t_total_full = t_cruise_full + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
L_D_cruise = 16.7
L_D_to = CL_to/CD_to
L_D_land = CL_land/ CD_to
tf =  0                     #Trap fuel time step
BSFC= 1/(43*10**6)   #Brake-specific fuel consumption
ddp = 0.8                   #Deep discharge protection
E_bat = 2.7*10**6           #Total Battery Energy per piece
eta_stt = 0.85 * 0.45       #Efficiency chain from shaft-to-thrust
eta_fuel_cell = 0.85 * 0.95 #Efficiency chain from shaft-to-thrust (fuel cell)
eta_btt = 0.95 * 0.85        #Efficiency chain from battery-to-thrust
NoD_ice = 2                 #Number of turboprop engines
NoD_em = 2                 #Number of electric motor engines

#Initial Climb
E_total_climb1 = (MTOW_design*72)/ (L_D_to) * t_climb1 + (MTOW_design/g * 72**2)/2 + MTOW_design*6.8*t_climb1
E_nc_climb1 = 0 * E_total_climb1     #Change %
E_c_climb1 = E_total_climb1 - E_nc_climb1
P_ice_climb1 = (E_c_climb1)/ (eta_stt * t_climb1 * NoD_ice)
m_fuel_climb1 = (1+tf)*P_ice_climb1*NoD_ice*BSFC*t_climb1
MTOW_climb_1 = MTOW_design - m_fuel_climb1 * g

#Second Climb
E_total_climb2 = (MTOW_climb_1*108)/ (L_D_cruise) * t_climb2 + (MTOW_climb_1/g * 36**2)/2 + MTOW_climb_1*5.08*t_climb2
E_nc_climb2 = 0 * E_total_climb2 #Change %
E_c_climb2 = E_total_climb2 - E_nc_climb2
P_ice_climb2 = (E_c_climb2)/ (eta_stt * t_climb2 * NoD_ice)
m_fuel_climb2 = (1+tf)*P_ice_climb2*NoD_ice*BSFC*t_climb2
MTOW_climb_2 = MTOW_climb_1 - m_fuel_climb2 * g

#3rd Climb
E_total_climb3 = (MTOW_climb_2*108)/ (L_D_cruise) * t_climb3 + MTOW_climb_2*4.06*t_climb3
E_nc_climb3 = 0 * E_total_climb3 #Change %
E_c_climb3 = E_total_climb3 - E_nc_climb3
P_ice_climb3 = (E_c_climb3)/ (eta_stt * t_climb3 * NoD_ice)
m_fuel_climb3 = (1+tf)*P_ice_climb3*NoD_ice*BSFC*t_climb3
MTOW_climb_3 = MTOW_climb_2 - m_fuel_climb3 * g

#Cruise
E_total_cruise_500 = (MTOW_climb_3*V_cruise)/ (L_D_cruise) * t_cruise_500 + (MTOW_climb_3/g * 33.4**2)/2
E_total_cruise_full = (MTOW_climb_3*V_cruise)/ (L_D_cruise) * t_cruise_full + (MTOW_climb_3/g * 33.4**2)/2
E_nc_cruise_500 = 0 * E_total_cruise_500 #Change %
E_c_cruise_500 = E_total_cruise_500 - E_nc_cruise_500
E_nc_cruise_full = 0 * E_total_cruise_full #Change %
E_c_cruise_full = E_total_cruise_full - E_nc_cruise_full
P_ice_cruise_500 = (E_c_cruise_500)/ (eta_stt * t_cruise_500 * NoD_ice)
P_ice_cruise_full = (E_c_cruise_full)/ (eta_stt * t_cruise_full * NoD_ice)
m_fuel_cruise_500 = (1+tf)*P_ice_cruise_500*NoD_ice*BSFC*t_cruise_500
m_fuel_cruise_full = (1+tf)*P_ice_cruise_full*NoD_ice*BSFC*t_cruise_full
MTOW_cruise_500 = MTOW_climb_3 - m_fuel_cruise_500 * g
MTOW_cruise_full = MTOW_climb_3 - m_fuel_cruise_full * g

#Descent 1
E_total_descent1 = (MTOW_cruise_full*138.9)/ (L_D_cruise) * t_descent1 + (MTOW_cruise_full/g * 25.7**2)/2 - MTOW_cruise_full*10.89*t_descent1
E_nc_descent1 = 0 * E_total_descent1 #Change %
E_c_descent1 = E_total_descent1 - E_nc_descent1
P_ice_descent1 = (E_c_descent1)/ (eta_stt * t_descent1 * NoD_ice)
m_fuel_descent1 = (1+tf)*P_ice_descent1*NoD_ice*BSFC*t_descent1
MTOW_descent1 = MTOW_cruise_full - m_fuel_descent1 * g

#Descent 2
E_total_descent2 = (MTOW_descent1*103)/ (L_D_land) * t_descent2 + (MTOW_descent1/g * 36**2)/2 - MTOW_descent1*8.08*t_descent2
E_nc_descent2 = 0 * E_total_descent2 #Change %
E_c_descent2 = E_total_descent2 - E_nc_descent2
P_ice_descent2 = (E_c_descent2)/ (eta_stt * t_descent2 * NoD_ice)
m_fuel_descent2 = (1+tf)*P_ice_descent2*NoD_ice*BSFC*t_descent2
MTOW_descent2 = MTOW_descent1 - m_fuel_descent2 * g

E_total_500 = E_total_climb1 + E_total_climb2 + E_total_climb3 + E_total_descent1 +E_total_descent2 +E_total_cruise_500
E_total_full = E_total_climb1 + E_total_climb2 + E_total_climb3 + E_total_descent1 +E_total_descent2 +E_total_cruise_full

# Hybridnization calculation
E_nc_total = E_nc_climb1 + E_nc_climb2 + E_nc_climb3 + E_nc_descent1 + E_nc_descent2 + E_nc_cruise_full
P_em = E_nc_total/ (eta_btt * t_total_full * NoD_em)
H_p_para = P_em*NoD_em / P_max
print(np.round(H_p_para,2)*100,"% Hybridlization")

# Mass Calculation
m_fuel_ice = m_fuel_climb1 + m_fuel_climb2 + m_fuel_climb3 + m_fuel_descent1 + m_fuel_descent2 + m_fuel_cruise_full
m_bat = (1+ddp) * (E_nc_total/(eta_btt*E_bat))
m_OE = (a * MTOW_design/g + b)
m_propulsion = m_turboprop*NoD_ice + (m_em + m_propeller)*NoD_em
m_MTOW = m_OE + m_fuel_ice + m_payload + m_bat + m_propulsion
print(np.round(m_fuel_ice,0), "Fuel Mass(kg)")
print(np.round(m_bat,0),"Battery Mass(kg)")
print(np.round(m_payload,0),"Payload Mass (kg)")
print(np.round(m_bat+m_fuel_ice,0),"Fuel + Battery mass")
print(np.round(m_OE,0), "Operational Empty mass without Engine(kg)")
print(np.round(m_propulsion,0), "Propulsion Mass (kg)")
print(np.round(m_MTOW,0), "Calculated MTOW (kg)")
print(np.round(MTOW_design/g,0), "MTOW Design(kg)")
print(np.round((m_MTOW - MTOW_design/g)/MTOW_design*100,2), "% Difference between hybrid and original design")
