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
WPAX = 200*0.453592*PAX*g                                #N
WPAXBAGGAGE = 40*0.453592*PAX*g                          #N Crew is bagageless
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

#Mass Preliminary Calculation
W_P_design = 0.0547
W_S_design = 3273
m_turboprop = 1074.5/2
m_em_dis = 13
m_propeller_dis = m_em_dis * 0.14
m_em_tip = 50
m_propeller_tip = m_em_tip * 0.14

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
MTOW_design = 20281 * g                  #N
S = MTOW_design / 3169
V_to = 1.13*(np.sqrt(1.1*MTOW_design/(1/2 * rho_1524 *S * CL_to)))
ROC1 = 4
ROC2 = 3
ROC3 = 2
ROD1 = 1500 / 196.9
ROD2 = 1110 / 196.9
kts_ms = 0.5144444
V_climb1 = 140 * (1 + 0.02 * 2.5) * kts_ms
V_climb2 = 176 * (1 + 0.02 * 10) * kts_ms#210 * kts_ms
V_climb3 = 176 * (1 + 0.02 * 21.5) * kts_ms#210 * kts_ms
V_descent1 = 176 * (1 + 0.02 * 19) * kts_ms#270 * kts_ms
V_descent2 = 140 * kts_ms
Vx_climb1 = np.sqrt(V_climb1**2 - ROC1**2)
Vx_climb2 = np.sqrt(V_climb2**2 - ROC2**2)
Vx_climb3 = np.sqrt(V_climb3**2 - ROC3**2)
Vx_descent1 = np.sqrt(V_descent1**2 - ROD1**2)
Vx_descent2 = np.sqrt(V_descent2**2 - ROD2**2)
t_climb1 = np.round(Climb1_h / ROC1,0)
a_climb1 = (V_climb1-V_to) / t_climb1
t_climb2 = np.round((Climb2_h - Climb1_h) / ROC2)
a_climb2 = (V_climb2 - V_climb1)/ t_climb2
t_climb3 = np.round((h_cruise - Climb2_h) / ROC3)
a_climb3 = (V_climb3 - V_climb2)/ t_climb3
t_descent1 = np.round((h_cruise - Descent1_h)/ ROD1)
t_descent2 = np.round((Descent1_h - Descent2_h)/ ROD2)
s_all_except_cruise= Vx_climb1 * t_climb1 + Vx_climb2 * t_climb2 + Vx_climb3 * t_climb3 + Vx_descent1 * t_descent1 + Vx_descent2 * t_descent2
s_cruise_500 = 500 * 1852 - s_all_except_cruise
s_cruise_full = R_eq - s_all_except_cruise
t_cruise_500 = np.ceil(s_cruise_500 / V_cruise)
t_cruise_full = np.ceil(s_cruise_full / V_cruise)
a_cruise_full = (V_cruise - V_climb3)/ t_cruise_full
a_cruise_500 = (V_cruise - V_climb3)/ t_cruise_500
a_descent1 = (V_descent1 - V_cruise) / t_descent1
a_descent2 = (V_descent2 - V_descent1) / t_descent2
t_total_500 = t_cruise_500 + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
t_total_full = t_cruise_full + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
L_D_cruise = CL/CD
L_D_to = CL_to/CD_to
L_D_land = CL_land/ CD_to
tf =  0                         #Trap fuel time step
BSFC= 1/(43*10**6 * 0.45)       #1/(43*10**6 * 0.39 * 0.9 *0.99)   #Brake-specific fuel consumption (only 43*10^6 * 0.45 if parallel series)
ddp = 0.8                       #Deep discharge protection
E_bat = 2.7*10**6               #Total Battery Energy per piece (Bat eff, Inverter eff, Em eff)
eta_stt = 0.85           #Efficiency chain from shaft-to-thrust
eta_btt = 0.934 * 0.85 * 0.99 * 0.995 * 0.95         #Efficiency chain from battery-to-thrust
NoD_ice = 2                     #Number of turboprop engines
NoD_em_tip = 2                  #Number of electric motor engines (wing tip)
NoD_em_dis = 2                  #Number of electric motor engines (distributed)
E_climb1_total = 0
E_climb2_total = 0
E_climb3_total = 0
E_cruise_full_total = 0
E_cruise_500_total = 0
E_descent1_total = 0
E_descent2_total = 0
E_nc_total = 0
m_fuel_total = 0
MTOW_energy_calculate = MTOW_design
hybridization_climb1 = []
hybridization_climb2 = []
hybridization_climb3 = []
hybridization_cruise_full = []
hybridization_descent1 = []
hybridization_descent2 = []
#Mass Preliminary Calculation
W_P_design = 0.0467
W_S_design = 3171
m_turboprop = 1074.5/2
m_em_dis = 13
m_propeller_dis = m_em_dis * 0.14
m_em_tip = 50
m_propeller_tip = m_em_tip * 0.14
P_max = MTOW_design / W_P_design
def m_fuel(P):
    m_fuel = P * BSFC
    return m_fuel
def P_ice(E,t):
    P_ice = E/(eta_stt * t)
    return P_ice
def Energy(MTOW,acc,t,L_D,ROC,V1):
    Energy = MTOW * (acc * t +V1) / L_D + (MTOW/g * acc **2)/2 + MTOW * ROC
    return Energy
def m_bat(E_nc):
    m_bat = (1+ddp) * (E_nc)/ (eta_btt * E_bat)
    return m_bat
def P_em(E_nc,dt):
    P_em = E_nc / (0.85 * dt)
    return P_em
for i in range(1,int(t_total_full)+1):
    if i<= t_climb1:
        E_climb1 = Energy(MTOW_energy_calculate,a_climb1,i,L_D_to,ROC1,V_to)
        E_nc = 0.55 * E_climb1
        E_c = E_climb1 - E_nc
        P_ice_climb1 = P_ice(E_c,1)
        P_em_climb1 = P_em(E_nc,1)
        P_max = P_ice_climb1 + P_em_climb1
        hybridization_climb1.append(P_em_climb1/P_max)
        m_fuel_climb1 = m_fuel(P_ice_climb1)
        E_nc_total += E_nc
        m_fuel_total += m_fuel_climb1
        E_climb1_total += E_climb1
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_climb1 * g
    elif t_climb1 < i <= (t_climb1+t_climb2):
        E_climb2 = Energy(MTOW_energy_calculate,a_climb2,(i-t_climb1),L_D_cruise,ROC2,V_climb1)
        E_nc = 0.55 * E_climb2
        E_c = E_climb2 - E_nc
        P_ice_climb2 = P_ice(E_c,1)
        P_em_climb2 = P_em(E_nc, 1)
        P_max = P_ice_climb2 + P_em_climb2
        hybridization_climb2.append(P_em_climb2/P_max)
        m_fuel_climb2 = m_fuel(P_ice_climb2)
        m_fuel_total += m_fuel_climb2
        E_climb2_total += E_climb2
        E_nc_total += E_nc
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_climb2 * g
    elif (t_climb1 + t_climb2) < i <= (t_climb1 +t_climb2 +t_climb3):
        E_climb3 = Energy(MTOW_energy_calculate,a_climb3,(i-(t_climb1+t_climb2)),L_D_cruise,ROC3,V_climb2)
        E_nc = 0.55 * E_climb3
        E_c = E_climb3 - E_nc
        P_ice_climb3 = P_ice(E_c, 1)
        P_em_climb3 = P_em(E_nc, 1)
        P_max = P_ice_climb3 + P_em_climb3
        hybridization_climb3.append(P_em_climb3 / P_max)
        m_fuel_climb3 = m_fuel(P_ice_climb3)
        m_fuel_total += m_fuel_climb3
        E_climb3_total += E_climb3
        E_nc_total += E_nc
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_climb3 * g
    elif (t_climb1 +t_climb2 +t_climb3) < i <= (t_climb1 + t_climb2 +t_climb3 + t_cruise_full):
        E_cruise_full = Energy(MTOW_energy_calculate, 0, (i -(t_climb1 +t_climb2 +t_climb3)), L_D_cruise, 0, V_cruise)
        E_nc =  0.0385 * E_cruise_full
        E_c = E_cruise_full - E_nc
        P_ice_cruise_full = P_ice(E_c, 1)
        P_em_cruise = P_em(E_nc, 1)
        P_max = P_ice_cruise_full + P_em_cruise
        hybridization_cruise_full.append(P_em_cruise / P_max)
        m_fuel_cruise = m_fuel(P_ice_cruise_full)
        E_nc_total += E_nc
        m_fuel_total += m_fuel_cruise
        E_cruise_full_total += E_cruise_full
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_cruise * g
    elif (t_climb1 + t_climb2 +t_climb3 + t_cruise_full) < i <= (t_climb1 + t_climb2 +t_climb3 + t_cruise_full + t_descent1):
        E_descent1 = Energy(MTOW_energy_calculate, a_descent1, (i - (t_climb1 + t_climb2 + t_climb3+t_cruise_full)), L_D_cruise,-ROD1, V_cruise)
        E_nc = 0 * E_descent1
        E_c = E_descent1 - E_nc
        P_ice_descent1 = P_ice(E_c, 1)
        P_em_descent1 = P_em(E_nc, 1)
        P_max = P_ice_descent1 + P_em_descent1
        hybridization_descent1.append(P_em_descent1 / P_max)
        m_fuel_descent1 = m_fuel(P_ice_descent1)
        m_fuel_total += m_fuel_descent1
        E_descent1_total += E_descent1
        E_nc_total += E_nc
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_descent1 * g
    else:
        E_descent2 = Energy(MTOW_energy_calculate, a_descent2, (i - (t_climb1 + t_climb2 +t_climb3 + t_cruise_full + t_descent1)), L_D_cruise, -ROD2, V_descent1)
        E_nc = 0 * E_descent2
        E_c = E_descent2 - E_nc
        P_ice_descent2 = P_ice(E_c, 1)
        P_em_descent2 = P_em(E_nc, 1)
        P_max = P_ice_descent2 + P_em_descent2
        hybridization_descent2.append(P_em_descent2 / P_max)
        m_fuel_descent2 = m_fuel(P_ice_descent2)
        m_fuel_total += m_fuel_descent2
        E_descent2_total += E_descent2
        E_nc_total += E_nc
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_descent2 * g

print(np.round(np.average(hybridization_climb1),2),"climb1 ")
print(np.round(np.average(hybridization_climb2),2),"climb2")
print(np.round(np.average(hybridization_climb3),2),"climb3")
print(np.round(np.average(hybridization_cruise_full),2),"Cruise")
print(np.round(np.average(hybridization_descent1),2),"descent1")
print(np.round(np.average(hybridization_descent2),2),"descent2")
m_propulsion = 700 * 1.95
m_OE = (a * MTOW_design/g + b) + m_propulsion
m_OE_without = (a * MTOW_design/g + b)
m_MTOW = m_OE + m_payload +m_bat(E_nc_total) + m_fuel_total
print(m_MTOW,"MTOM")
print(m_OE, "OEM")
print(m_bat(E_nc_total),"m_bat")
print(m_fuel_total,"m_fuel")
S = m_MTOW*g / W_S_design
print(S,"m^2")
b = np.sqrt(A*S)
print(b,"m")
print(m_payload)
print(m_OE_without)
print((E_climb1_total + E_climb2_total + E_climb3_total + E_cruise_full_total + E_descent1_total + E_descent2_total)/10**6)