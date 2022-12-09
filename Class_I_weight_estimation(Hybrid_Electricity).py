import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
lambda_trop = -6.5/1000
R = 287.0528
V_cruise = 275 * 0.51444444     #knt -> m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
s_takeoff_1524 = 1372           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 1372           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
rho_0= 1.225                    #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524= 1.01893               #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
eff_prop = 0.85                 #Change with Literature
PAX = 50
WPAX = 200*0.453592*PAX*g + 3*190*0.453592*g                                  #N
WPAXBAGGAGE = 40*0.453592*PAX*g +2*30*0.453592*g                              #N Crew is bagageless
m_payload = (WPAX + WPAXBAGGAGE) / g
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
TOP = 430                          #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC = 6.9                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0032                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach = 141* 0.514444         #Change with CS25 or Requirement
a = 0.5088
b = 1199.7

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to #* rho_1524_rho0 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_1524 * s_landing_1524/0.5847)/(2*0.95)
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_cruise_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + (CD_to/CL_to))*(np.sqrt((2/rho_0)*(1/CL_to))))
#Stall Constraint
W_S_approach = 1/2 * rho_1524 * V_approach**2 * CL_land

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
plt.title("Loading Diagram of Novel Regional Aircraft")
plt.legend(loc = "upper right")
plt.grid()
plt.show()

#Mass Preliminary Calculation
W_P_design = 0.0547
W_S_design = 3273
m_turboprop = 1074.5/2
m_em = 13
m_propeller = m_em * 0.14
MTOW_design = 21123 * g                  #N

#Range Calculation
CL = np.sqrt(np.pi*Cd0*A*e)
CD = CD = 2 * Cd0
R_norm = 1000 * 1852
R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))
f_con = 0.05
R_div = 100 * 1852
E = 45 * 60
R_eq = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + (E*V_cruise)
#Energy Calculation
Climb1_h = 50 * 100 *0.3048
Climb2_h = 150 * 100 *0.3048
Descent1_h = 100 * 100 *0.3048
Descent2_h = 0
ROC1 = 1350 / 196.9
ROC2 = 1000 / 196.9
ROC3 = 800 / 196.9
ROD1 = 1500 / 196.9
ROD2 = 1110 / 196.9
V_climb1 = 140 * 0.51444444
V_climb2 = 210 * 0.51444444
V_climb3 = 210 * 0.51444444
V_descent1 = 270 * 0.51444444
V_descent2 = 140 * 0.51444444
Vx_climb1 = np.sqrt(V_climb1**2 - ROC1**2)
Vx_climb2 = np.sqrt(V_climb2**2 - ROC2**2)
Vx_climb3 = np.sqrt(V_climb3**2 - ROC3**2)
Vx_descent1 = np.sqrt(V_descent1**2 - ROD1**2)
Vx_descent2 = np.sqrt(V_descent2**2 - ROD2**2)
t_climb1 = Climb1_h / ROC1
t_climb2 = (Climb2_h - Climb1_h) / ROC2
t_climb3 = (h_cruise - Climb2_h) / ROC3
t_descent1 = (h_cruise - Descent1_h)/ ROD1
t_descent2 = (Descent1_h - Descent2_h)/ ROD2
s_all_except_cruise= Vx_climb1 * t_climb1 + Vx_climb2 * t_climb2 + Vx_climb3 * t_climb3 + Vx_descent1 * t_descent1 + Vx_descent2 * t_descent2
s_cruise_500 = 500 * 1852 - s_all_except_cruise
s_cruise_full = R_eq - s_all_except_cruise
t_cruise_500 = s_cruise_500 / V_cruise
t_cruise_full = s_cruise_full / V_cruise
t_total_500 = t_cruise_500 + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
t_total_full = t_cruise_full + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
L_D_cruise = CL/CD
L_D_to = CL_to/CD_to
L_D_land = CL_land/ CD_to
tf =  0                     #Trap fuel time step
BSFC= 1/(43*10**6 * 0.45)   #Brake-specific fuel consumption
ddp = 0.8                   #Deep discharge protection
E_bat = 2.7*10**6           #Total Battery Energy per piece
eta_stt = 0.85       #Efficiency chain from shaft-to-thrust
eta_btt = 0.934 * 0.85        #Efficiency chain from battery-to-thrust
NoD_ice = 2                   #Number of turboprop engines
NoD_em = 2                 #Number of electric motor engines

#Initial Climb
E_total_climb1 = (MTOW_design*V_climb1)/ (L_D_to) * t_climb1 + (MTOW_design/g * V_climb1**2)/2 + MTOW_design*ROC1*t_climb1
E_nc_climb1 = 0 * E_total_climb1     #Change %
E_c_climb1 = E_total_climb1 - E_nc_climb1
P_ice_climb1 = (E_c_climb1)/ (eta_stt * t_climb1 * NoD_ice)
m_fuel_climb1 = (1+tf)*P_ice_climb1*NoD_ice*BSFC*t_climb1
MTOW_climb_1 = MTOW_design - m_fuel_climb1 * g

#Second Climb
E_total_climb2 = (MTOW_climb_1*V_climb2)/ (L_D_cruise) * t_climb2 + (MTOW_climb_1/g * (V_climb2-V_climb1)**2)/2 + MTOW_climb_1*ROC2*t_climb2
E_nc_climb2 = 0 * E_total_climb2 #Change %
E_c_climb2 = E_total_climb2 - E_nc_climb2
P_ice_climb2 = (E_c_climb2)/ (eta_stt * t_climb2 * NoD_ice)
m_fuel_climb2 = (1+tf)*P_ice_climb2*NoD_ice*BSFC*t_climb2
MTOW_climb_2 = MTOW_climb_1 - m_fuel_climb2 * g

#3rd Climb
E_total_climb3 = (MTOW_climb_2*V_climb3)/ (L_D_cruise) * t_climb3 + MTOW_climb_2*ROC3*t_climb3
E_nc_climb3 = 0 * E_total_climb3 #Change %
E_c_climb3 = E_total_climb3 - E_nc_climb3
P_ice_climb3 = (E_c_climb3)/ (eta_stt * t_climb3 * NoD_ice)
m_fuel_climb3 = (1+tf)*P_ice_climb3*NoD_ice*BSFC*t_climb3
MTOW_climb_3 = MTOW_climb_2 - m_fuel_climb3 * g

#Cruise
E_total_cruise_500 = (MTOW_climb_3*V_cruise)/ (L_D_cruise) * t_cruise_500 + (MTOW_climb_3/g * (V_cruise-V_climb3)**2)/2
E_total_cruise_full = (MTOW_climb_3*V_cruise)/ (L_D_cruise) * t_cruise_full + (MTOW_climb_3/g * (V_cruise-V_climb3)**2)/2
E_nc_cruise_500 = 0 * E_total_cruise_500 #Change %
E_c_cruise_500 = E_total_cruise_500 - E_nc_cruise_500
E_nc_cruise_full = E_nc_cruise_500 #0 * E_total_cruise_full #Change %
E_c_cruise_full = E_total_cruise_full - E_nc_cruise_500
P_ice_cruise_500 = (E_c_cruise_500)/ (eta_stt * t_cruise_500 * NoD_ice)
P_ice_cruise_full = (E_c_cruise_full)/ (eta_stt * t_cruise_full * NoD_ice)
m_fuel_cruise_500 = (1+tf)*P_ice_cruise_500*NoD_ice*BSFC*t_cruise_500
m_fuel_cruise_full = (1+tf)*P_ice_cruise_full*NoD_ice*BSFC*t_cruise_full
MTOW_cruise_500 = MTOW_climb_3 - m_fuel_cruise_500 * g
MTOW_cruise_full = MTOW_climb_3 - m_fuel_cruise_full * g

#Descent 1
E_total_descent1 = (MTOW_cruise_full*V_descent1)/ (L_D_cruise) * t_descent1 + (MTOW_cruise_full/g * (V_cruise-V_descent1)**2)/2 - MTOW_cruise_full*ROD1*t_descent1
E_nc_descent1 = 0 * E_total_descent1 #Change %
E_c_descent1 = E_total_descent1 - E_nc_descent1
P_ice_descent1 = (E_c_descent1)/ (eta_stt * t_descent1 * NoD_ice)
m_fuel_descent1 = (1+tf)*P_ice_descent1*NoD_ice*BSFC*t_descent1
MTOW_descent1 = MTOW_cruise_500 - m_fuel_descent1 * g

#Descent 2
E_total_descent2 = (MTOW_descent1*V_descent2)/ (L_D_land) * t_descent2 + (MTOW_descent1/g * (V_descent1 - V_descent2)**2)/2 - MTOW_descent1*ROD2*t_descent2
E_nc_descent2 = 0 * E_total_descent2 #Change %
E_c_descent2 = E_total_descent2 - E_nc_descent2
P_ice_descent2 = (E_c_descent2)/ (eta_stt * t_descent2 * NoD_ice)
m_fuel_descent2 = (1+tf)*P_ice_descent2*NoD_ice*BSFC*t_descent2
MTOW_descent2 = MTOW_descent1 - m_fuel_descent2 * g

E_total_500 = E_total_climb1 + E_total_climb2 + E_total_climb3 + E_total_descent1 +E_total_descent2 +E_total_cruise_500
E_total_full = E_total_climb1 + E_total_climb2 + E_total_climb3 + E_total_descent1 +E_total_descent2 +E_total_cruise_full

# Hybridnization calculation
E_nc_total = E_nc_climb1 + E_nc_climb2 + E_nc_climb3 + E_nc_descent1 + E_nc_descent2 + E_nc_cruise_500
P_em = E_nc_total/ (eta_btt * t_total_full * NoD_em)

#Mass Calculation
m_fuel_ice = m_fuel_climb1 + m_fuel_climb2 + m_fuel_climb3 + m_fuel_descent1 + m_fuel_descent2 + m_fuel_cruise_full
m_bat = (1+ddp) * (E_nc_total/(0.95*E_bat))
m_propulsion = m_turboprop*NoD_ice + (m_em + m_propeller)*NoD_em
m_OE = (a * MTOW_design/g + b) + m_propulsion
m_OE_without = (a * MTOW_design/g + b)
m_MTOW = m_OE + m_fuel_ice + m_payload + m_bat
print(np.round(m_fuel_ice,0), "Fuel Mass(kg)")
print(np.round(m_bat,0),"Battery Mass(kg)")
print(np.round(m_payload,0),"Payload Mass (kg)")
print(np.round(m_bat+m_fuel_ice,0),"Fuel + Battery mass")
print(np.round(m_OE,0), "Operational Empty mass with Engine(kg)")
print(np.round(m_OE_without,0), "Operational Empty mass without Engine(kg)")
print(np.round(m_propulsion,0), "Propulsion Mass (kg)")
print(np.round(m_MTOW,0), "Calculated MTOM (kg)")
print(np.round(MTOW_design/g,0), "Initial MTOM Input(kg)")
print(np.round((m_MTOW/(MTOW_design/g))-1,2), "Difference between hybrid and original design")
P_max = m_MTOW*g / W_P_design
S = m_MTOW*g / W_S_design
print(P_max/1000, "kW Max Power")
print(S,"m^2 Surface Area ")
H_p_para = P_em*NoD_em / P_max
print(np.round(H_p_para,2)*100,"% Hybridlization")