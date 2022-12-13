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
WPAX = 200*0.453592*PAX*g                                 #N
WPAXBAGGAGE = 40*0.453592*PAX*g                           #N Crew is bagageless
m_payload = (WPAX + WPAXBAGGAGE) / g
W_S = np.arange(1,4000,1)

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
ROC = 4.07                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach = 141* 0.514444         #Change with CS25 or Requirement
a = 0.5088
b = 1199.7

#Mass Preliminary Calculation
W_P_design_cruise = 0.04666
W_P_design_roc1 = 0.07426
W_P_design_roc2 = 0.08830
W_P_design_roc3 = 0.09867
W_S_design = 3273
m_turboprop = 1074.5/2
m_em_dis = 13
m_propeller_dis = m_em_dis * 0.14
m_em_tip = 50
m_propeller_tip = m_em_tip * 0.14
MTOW_design = 21123 * g                  #N
S = MTOW_design/ W_S_design
#Range Calculation
CL = np.sqrt(np.pi*Cd0*A*e)
CD = 2 * Cd0
R_norm = 1000 * 1852
R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))
f_con = 0.05
R_div = 100 * 1852
E = 30 * 60
R_eq = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + (E*V_cruise)

P_cruise_max = MTOW_design / W_P_design_cruise
P_climb1_max = MTOW_design / W_P_design_roc1
P_climb2_max = MTOW_design / W_P_design_roc2
P_climb3_max = MTOW_design / W_P_design_roc3
'''print(P_cruise_max/1000)
print(P_climb1_max/1000)
print(P_climb2_max/1000)
print(P_climb3_max/1000)'''
#Energy Calculation
Climb1_h = 50 * 100 *0.3048
Climb2_h = 150 * 100 *0.3048
Descent1_h = 100 * 100 *0.3048
Descent2_h = 0
V_to = 1.13*(np.sqrt(1.1*MTOW_design/(1/2 * rho_1524 *S * CL_to)))
ROC1 = 1350 / 196.9
ROC2 = 1000 / 196.9
ROC3 = 800 / 196.9
ROD1 = 1500 / 196.9
ROD2 = 1110 / 196.9
V_climb1 = 140 * 0.51444444
V_climb2 = 210 * 0.51444444
V_climb3 = 210 * 0.51444444
V_descent1 = 270 * 0.51444444
V_descent2 = 141 * 0.51444444
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
L_D_land = CL_land/ CD_land
tf =  0                     #Trap fuel time step
BSFC= 1/(43*10**6 * 0.39 * 0.9 *0.99)   #Brake-specific fuel consumption (only 43*10^6 * 0.45 if parallel series)
ddp = 0.8                   #Deep discharge protection
E_bat = 2.7*10**6      #Total Battery Energy per piece
eta_stt = 0.85                #Efficiency chain from shaft-to-thrust
eta_btt = 0.934 * 0.85 *0.99 * 0.995 * 0.95        #Efficiency chain from battery-to-thrust (Bat eff, Inverter eff, Em eff)
NoD_ice = 2                   #Number of turboprop engines
NoD_em_tip = 2                #Number of electric motor engines (wing tip)
NoD_em_dis = 0                #Number of electric motor engines (distributed)
def m_fuel(P_ice,dt):
    m_fuel = (1+tf)*P_ice*NoD_ice*BSFC*dt
    return m_fuel
def E_nc(P_em,dt):
    E_nc = P_em * dt * (NoD_em_tip) * eta_stt
    return E_nc
def P_em(P_max,doh):
    P_em = P_max * doh
    return P_em
def m_bat(E_nc):
    m_bat = (1+ddp) * (E_nc)/ (eta_btt * E_bat)
    return m_bat
#Takeoff


#Initial Climb
doh_climb1 = 0.55
E_climb_total = P_cruise_max * (t_climb1 * eta_stt)
print(E_climb_total/10**6)
P_em_climb1 = P_em(P_cruise_max,doh_climb1)
P_ice_climb1 = P_cruise_max - P_em_climb1
E_nc_climb1 = E_nc(P_em_climb1,t_climb1)
m_fuel_climb1 = m_fuel(P_ice_climb1,t_climb1)
m_bat_climb1 = m_bat(E_nc_climb1)

#Climb2
doh_climb2 = 0.55
P_em_climb2 = P_em(P_cruise_max,doh_climb2)
E_nc_climb2 = E_nc(P_em_climb2,t_climb2)
P_ice_climb2 = P_cruise_max - P_em_climb1
m_fuel_climb2 = m_fuel(P_ice_climb2,t_climb2)
m_bat_climb2 = m_bat(E_nc_climb2)

#Climb3
doh_climb3 = 0.55
P_em_climb3 = P_em(P_cruise_max,doh_climb3)
P_ice_climb3 = P_cruise_max - P_em_climb3
E_nc_climb3 = E_nc(P_em_climb3,t_climb3)
m_fuel_climb3 = m_fuel(P_ice_climb3,t_climb3)
m_bat_climb3 = m_bat(E_nc_climb3)

#Cruise
doh_cruise = 0.2
P_em_cruise = P_em(P_cruise_max,doh_cruise)
P_ice_cruise = P_cruise_max - P_em_cruise
E_nc_cruise = E_nc(P_em_cruise,t_cruise_full)
m_fuel_cruise = m_fuel(P_ice_cruise,t_cruise_full)
m_bat_cruise = m_bat(E_nc_cruise)

#Descent1
doh_descent1 = 0
P_em_cruise = P_em(P_cruise_max,doh_descent1)

#Descent2
doh_descent2 = 0
P_em_cruise = P_em(P_cruise_max,doh_descent2)
print(m_bat_climb1)
print(m_fuel_climb1)
print(m_bat_climb2)
print(m_fuel_climb2)
print(m_bat_climb3)
print(m_fuel_climb3)
print(m_bat_cruise)
print(m_fuel_cruise)
