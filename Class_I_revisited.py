import numpy as np
import matplotlib.pyplot as plt

#Constants
g = 9.80665
Molar_mass_air = 0.0289644 #kg/mol
universal_gas_constant = 8.31432 #N m kmol⁻¹ K⁻¹
specific_gas_constant = 287.052 #J·kg⁻¹·K⁻¹
kts_ms = 0.514444444
ft_to_m = 0.3048

T_0 = 288.15
p_0 = 101325 #Pa
rho_0 = 1.225
lapse_rate = -0.0065

V_cruise = 275 * kts_ms #KTAS
h_cruise = 28000 * ft_to_m

s_takeoff = 1372 #m
s_landing = 1372 #m
alt_1524 = 5000 * ft_to_m

def ISA_calculator(h):
    T = T_0 + lapse_rate * h
    p = p_0 * ((T_0 / T) ** ((g * Molar_mass_air) / (universal_gas_constant * lapse_rate)))
    rho = p / (specific_gas_constant * T)
    return T, p, rho

T_1524, p_1524, rho_1524 = ISA_calculator(alt_1524)
T_cruise, p_cruise, rho_cruise = ISA_calculator(h_cruise)
sigma_1524 = rho_1524 / rho_0
sigma_cruise = rho_cruise / rho_0

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
ROC = 6.9                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach = 141* kts_ms         #Change with CS25 or Requirement
a = 0.5088
b = 1199.7

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to * sigma_1524 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_1524 * s_landing/0.5847)/(2*0.95)
# Cruise Speed Constraint
W_P_cru = eff_prop * (sigma_cruise)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + (CD_to/CL_to))*(np.sqrt((2/rho_0)*(1/CL_to))))
#Stall Constraint
W_S_approach = 1/2 * rho_1524 * V_approach**2 * CL_land


plt.vlines(W_S_approach,0,100,'b',label="Approach Speed Constraint")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,100,'k',label ="Landing Constraint")
plt.axvspan(3271,4000,color = "red", alpha = 0.1)
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.fill_between(W_S,W_P_cru,1,color = "red",alpha = 0.1)
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,4000)
plt.ylim(0,0.5)
plt.xticks(np.arange(0,4001,500))
plt.yticks(np.arange(0,0.5,0.05))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.legend(loc = "upper right")
plt.grid()
plt.show()



#Mass Preliminary Calculation
W_P_design = 0.04706
W_S_design = 3271
m_turboprop = 1074.5/2
MTOW_design = 20281 * g                  #N

#Range Calculation
CL = np.sqrt(np.pi*Cd0*A*e)
CD = 2 * Cd0
R_norm = 1000 * 1852
R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))
f_con = 0.05
R_div = 100 * 1852
E = 45 * 60
R_eq = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + (E*V_cruise)

climb_h1 = 5000 * ft_to_m
climb_h2 = 15000 * ft_to_m - climb_h1
climb_h3 = 28000 * ft_to_m - climb_h2 - climb_h1
descent_h1 = climb_h3 - 10000 * ft_to_m
descent_h2 = 10000 * ft_to_m
h_all = [climb_h1, climb_h2, climb_h3, descent_h1, descent_h2]

v_climb1 = 140 * (1 + 0.02 * 2.5) * kts_ms
v_climb2 = 176 * (1 + 0.02 * 10) * kts_ms#210 * kts_ms
v_climb3 = 176 * (1 + 0.02 * 21.5) * kts_ms#210 * kts_ms
v_descent1 = 176 * (1 + 0.02 * 19) * kts_ms#270 * kts_ms
v_descent2 = 140 * kts_ms
V_all = [v_climb1, v_climb2, v_climb3, v_descent1, v_descent2]
print(V_all)
ROC1 = 4#1350 / 196.9
ROC2 = 3#1000 / 196.9
ROC3 = 2#800 / 196.9
ROD1 = 1500 / 196.9
ROD2 = 1110 / 196.9
ROC_all = [ROC1, ROC2, ROC3, ROD1, ROD2]

def t_cruise():
    V_x = []
    t_list = []
    for i in np.arange(0, len(ROC_all)):
        V_x.append(np.sqrt(V_all[i]**2 - ROC_all[i]**2))
        t_list.append(h_all[i] / ROC_all[i])
    s_no_cruise = np.array(V_x) * np.array(t_list)
    s_cruise_500 = 500 * 1852 - sum(s_no_cruise)
    s_cruise_full = R_eq - sum(s_no_cruise)
    t_cruise_500 = s_cruise_500 / V_cruise
    t_cruise_full = s_cruise_full / V_cruise
    return t_cruise_500, t_cruise_full, t_list

t_cruise_500, t_cruise_full, t_list = t_cruise()
t_total_500 = sum(t_list) + t_cruise_500 
t_total_full = sum(t_list) + t_cruise_full

L_D_cruise = CL/CD
L_D_to = CL_to/CD_to
L_D_land = CL_land/ CD_to
tf =  0                     #Trap fuel time step
BSFC= 1/(43*10**6)   #Brake-specific fuel consumption
ddp = 0.8                   #Deep discharge protection
BSFC_lh2 = 1/ (120*10**6)           #Total Battery Energy per piece
eta_stt = 0.39 * 0.9 *0.99                 #Efficiency chain from shaft-to-thrust
eta_fuel_cell = 0.6 * 0.97 * 0.995**2 * 0.85 * 0.95    #Efficiency chain from shaft-to-thrust (fuel cell)
NoD_ice = 2                 #Number of turboprop engines

#Power Calculation
V_to = 65 #m/s
E_total = []
P_fc_total = []
MTOW_total = []
m_fc_total = []
def power(V, V_prev, t, ROC, MTOW, L_D):
    E_increment = (MTOW * V * t) / (L_D) + (MTOW / g * (V - V_prev) ** 2) / 2 + MTOW * ROC * t
    P_fc = E_increment / (eta_fuel_cell * t * NoD_ice)
    m_fc = (1+tf) * P_fc * NoD_ice * BSFC_lh2 * t
    MTOW_now = MTOW - m_fc * g
    E_total.append(E_increment / 10**6)
    P_fc_total.append(P_fc / 10**3)
    MTOW_total.append(MTOW_now)
    m_fc_total.append(m_fc)
    return E_total, P_fc_total, MTOW_total, m_fc_total

def total_all(range):
    power(v_climb1, V_to, t_list[0], ROC1, MTOW_design, L_D_to)
    power(v_climb2, v_climb1, t_list[1], ROC2, MTOW_total[0], L_D_cruise)
    power(v_climb3, v_climb2, t_list[2], ROC3, MTOW_total[1], L_D_cruise)
    if range == 500:
        power(V_cruise, v_climb3, t_cruise_500, 0, MTOW_total[2], L_D_cruise)
    else:
        power(V_cruise, v_climb3, t_cruise_full, 0, MTOW_total[2], L_D_cruise)
    power(v_descent1, V_cruise, t_list[3], -1 * ROD1, MTOW_total[3], L_D_cruise)
    power(v_descent2, v_descent1, t_list[4], -1 * ROD2, MTOW_total[4], L_D_land)
    return E_total, P_fc_total, MTOW_total, m_fc_total

def sum(range):
    total_all(range)
    E_sum = np.sum(E_total)
    P_sum = np.sum(P_fc_total)
    MTOW_sum = np.sum(MTOW_total)
    m_fc_sum = np.sum(m_fc_total)
    return E_sum, P_sum, MTOW_sum, m_fc_sum
# E_sum_500, P_sum_500, MTOW_sum_500 = sum(500)
E_sum_full, P_sum_full, MTOW_sum_full, m_fc_sum_full = sum(1000)
P_max = max(P_fc_total) #in kW
m_fuelcell_struc = (P_max / 1000)/3
m_inverter = (P_max)/30
m_propulsion = (m_turboprop*NoD_ice + m_inverter)* 1.5 + m_fuelcell_struc
m_propulsion_withoutstruc = m_propulsion - m_fuelcell_struc
m_OE = (a * MTOW_design/g + b) + m_propulsion
m_OE_without = (a * MTOW_design/g + b)
m_MTOW = m_OE + m_fc_sum_full + m_payload

print(f"Total energy is {E_sum_full} MJ")
print(f"Peak power is {P_max} kW per engine")
print(P_fc_total, t_list, t_cruise_full)
print(f"The MTOW is {m_MTOW} kg")
print(E_total)
print(m_fc_total)




