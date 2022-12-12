import numpy as np
#Energy Calculation
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
eff_prop = 0.85                 #Change with Literature

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
a_climb1 = V_climb1 / t_climb1
t_climb2 = (Climb2_h - Climb1_h) / ROC2
a_climb2 = (V_climb2 - V_climb3)/ t_climb1
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
MTOW_design = 20281 * g                  #N
E_climb1 = 0
print(t_total_full)
eta_stt = 0.85
for s in range(1,223):
    if s <= t_climb1:
        E_climb1 += (MTOW_design*(a_climb1*s))/ (L_D_to) + (MTOW_design/g * a_climb1**2)/2 + MTOW_design*ROC1
    else:
        continue
P_climb1 = E_climb1/(0.45 * t_climb1)
a_cruise = (V_cruise - V_climb3)/t_cruise_full
E_cruise = 0
MTOW_climb3 = 198163.79275823693
for s in range(1,17358):
    if s<= t_cruise_full:
        E_cruise += (MTOW_climb3*(a_cruise*s + V_climb3))/ (L_D_cruise) + (MTOW_climb3/g * a_cruise**2)/2
    else:
        continue
P_cruise = E_cruise/(0.45 *t_cruise_full)
print(E_climb1/10**6,"Energy Climb1")
print(P_climb1/1000, "Power Climb1")
print(E_cruise/10**6, "Energy Cruise")
print(P_cruise/1000, "Power Cruise")