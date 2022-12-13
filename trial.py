import numpy as np
#Constants
g = 9.80665
lambda_trop = -6.5/1000
R = 287.0528
V_cruise = 275 * 0.51444444     #knt -> m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
rho_0= 1.225                    #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524= 1.01893               #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
##Cdo calculations
Psi = 0.0075                    #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #span efficiency factor (value based on Roelof reader p.46)
A = 12                           #Aspect Ratio (12-14) #Reference to ATR 72
e = 1/(np.pi*A*Psi+(1/phi))
Cfe = 0.0045                     #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1                     #(6.0-6.2) wetted area ratios -> depending on airframe structure
Cd0 = Cfe * Swet_S
#Range Calculation
CL = np.sqrt(np.pi*Cd0*A*e)
CD = CD = 2 * Cd0
R_norm = 1000 * 1852
R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))
f_con = 0.05
R_div = 100 * 1852
E = 45 * 60
R_eq = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + (E*V_cruise)
#Aerodynamic Estimations
CL_max = 1.9                       #(1.2-1.8 twin engine) & (1.5-1.9 turboprop) max lift coefficient
CL_to = 2.1                        #Change with Estimate (1.7-2.1)
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
CL_land = 2.6                      #Change with Estimate (1.9-3.3)
CD_land = Cd0 + (CL_to**2 /(np.pi * A* e))
TOP = 430                          #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC = 6.9                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
Climb1_h = 50 * 100 *0.3048
Climb2_h = 150 * 100 *0.3048
Descent1_h = 100 * 100 *0.3048
Descent2_h = 0
MTOW_design = 20281 * g                  #N
S = MTOW_design / 3169
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
tf =  0                     #Trap fuel time step
BSFC= 1/(43*10**6 * 0.45) #1/(43*10**6 * 0.39 * 0.9 *0.99)   #Brake-specific fuel consumption (only 43*10^6 * 0.45 if parallel series)
ddp = 0.8                   #Deep discharge protection
E_bat = 2.7*10**6 *0.99 * 0.995 * 0.95      #Total Battery Energy per piece (Bat eff, Inverter eff, Em eff)
eta_stt = 0.85                #Efficiency chain from shaft-to-thrust
eta_btt = 0.934 * 0.85        #Efficiency chain from battery-to-thrust
NoD_ice = 2                   #Number of turboprop engines
NoD_em_tip = 2                #Number of electric motor engines (wing tip)
NoD_em_dis = 2                #Number of electric motor engines (distributed)
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
def m_fuel(P):
    BSFC = 1/ (43 * 10**6 * 0.45)
    m_fuel = P * BSFC
    return m_fuel
def P_ice(E,t):
    P_ice = E/(eta_stt * t)
    return P_ice
def Energy(MTOW,acc,t,L_D,ROC,V1):
    Energy = MTOW * (acc * t +V1) / L_D + (MTOW/g * acc **2)/2 + MTOW * ROC
    return Energy
for i in range(1,int(t_total_full)+1):
    if i<= t_climb1:
        E_climb1 = Energy(MTOW_energy_calculate,a_climb1,i,L_D_to,ROC1,V_to)
        E_nc = 0 * E_climb1
        E_c = E_climb1 - E_nc
        P_ice_climb1 = P_ice(E_c,1)
        m_fuel_climb1 = m_fuel(P_ice_climb1)
        E_nc_total = E_nc
        m_fuel_total = m_fuel_climb1
        E_climb1_total += E_climb1
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_climb1 * g
    elif t_climb1 < i <= (t_climb1+t_climb2):
        E_climb2 = Energy(MTOW_energy_calculate,a_climb2,(i-t_climb1),L_D_cruise,ROC2,V_climb1)
        E_nc = 0 * E_climb2
        E_c = E_climb2 - E_nc
        P_ice_climb2 = P_ice(E_c,1)
        m_fuel_climb2 = m_fuel(P_ice_climb2)
        m_fuel_total = m_fuel_climb2
        E_climb2_total += E_climb2
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_climb2 * g
    elif (t_climb1 + t_climb2) < i <= (t_climb1 +t_climb2 +t_climb3):
        E_climb3 = Energy(MTOW_energy_calculate,a_climb3,(i-(t_climb1+t_climb2)),L_D_cruise,ROC3,V_climb2)
        E_nc = 0 * E_climb3
        E_c = E_climb3 - E_nc
        P_ice_climb3 = P_ice(E_c, 1)
        m_fuel_climb3 = m_fuel(P_ice_climb3)
        m_fuel_total = m_fuel_climb3
        E_climb3_total += E_climb3
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_climb3 * g
    elif (t_climb1 +t_climb2 +t_climb3) < i <= (t_climb1 + t_climb2 +t_climb3 + t_cruise_full):
        E_cruise_full = Energy(MTOW_energy_calculate,a_cruise_full,(i-(t_climb1+t_climb2+t_climb3)),L_D_cruise,0,V_climb3)
        E_nc = 0 * E_cruise_full
        E_c = E_cruise_full - E_nc
        P_ice_cruise_full = P_ice(E_c,1)
        m_fuel_cruise_full = m_fuel(P_ice_cruise_full)
        E_cruise_full_total += E_cruise_full
        MTOW_energy_calculate = MTOW_energy_calculate - m_fuel_cruise_full *g
        continue
print(E_climb1_total/10**6,"MJ")
print(E_climb1_total/t_climb1/1000,"kW")
print(E_climb2_total/10**6, "MJ")
print(E_climb2_total/t_climb2/1000,"kW")
print(E_climb3_total/10**6, "MJ")
print(E_climb3_total/t_climb3/1000,"kW")
print(E_cruise_full_total/10**6, "MJ")
print(E_cruise_full_total/t_cruise_full/1000,"kW")
print(MTOW_energy_calculate)
print(MTOW_design)
print(t_cruise_full)
print(L_D_cruise)
print(h_cruise)
print(V_climb3)
print(a_cruise_full)