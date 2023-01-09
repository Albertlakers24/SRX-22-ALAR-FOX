import numpy as np
from matplotlib import pyplot as plt
from Constants import *
#Constants
rho_1524= ISA_calculator(4500*ft_m,dt_takeoff)[2]              #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise = ISA_calculator(280*FL_ft*ft_m,dt_cruise)[2]
W_S = np.arange(1,6000,1)
s_takeoff_1524 = s_takeoff
s_landing_1524 = s_landing
##Cdo calculations
e = 1/(np.pi*A*Psi+(1/phi))
Cd0 = Cfe * Swet_S
#ROC and beta estimates
ROC = 4                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.024#0.0032
ROC_V_OEI = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach_stall = V_approach /1.23  #CS 25 requirement of V_stall_land = V_approach / 1.23
beta_V_app_fc = 0.97
beta_s_land_fc = 0.97
beta_cruise_fc  = 0.985
beta_ROC_fc = 0.985
beta_V_app_hc = 0.94
beta_s_land_hc = 0.94
beta_cruise_hc  = 0.97
beta_ROC_hc = 0.97
beta_V_app_ker = 0.91
beta_s_land_ker = 0.91
beta_cruise_ker  = 0.96
beta_ROC_ker = 0.96
beta_cv = 1
beta_s_to = 1
beta_em = 1
C_LFL = 0.45                        #Landing field length coefficient s^2/m
alpha_p_em = 1
CL_2 = 1.13**2 * CL_max_takeoff
k_t = 0.85
h2 = 50 * ft_m
V_s0 = 31
def oswald_efficiency(flap_deflection):
    delta_e = 0.0026 * flap_deflection
    e_new = e + delta_e
    return e_new
def CD_0(flap_deflection, lg):
    delta_CD0 = flap_deflection * (13 * 10**(-4)) + (175 * 10**(-4)) * lg
    CD_0_new = Cd0 + delta_CD0
    return CD_0_new
def Power_lapse(rho,rho_0):
    alpha_p_ice = (rho/rho_0) ** (3/4)
    return alpha_p_ice
def V_approach_constraint(density,V_approach_stall,CL_max_landing,beta):
    W_S = 1/2 * density * V_approach_stall **2 * CL_max_landing * 1/beta
    return W_S
def s_land_constraint(s_land,C_LFL,density,CL_max_landing,beta):
    W_S = (s_land / C_LFL) * (density * CL_max_landing / 2) * 1/beta
    return W_S
def cruise_contraint(eta_prop,alpha_p,Cd0,density,V,W_S,A,e,beta):
    W_P = eta_prop * (alpha_p/beta) *((Cd0*1/2*density*V**3)/(beta*W_S)+ (beta*W_S)/(np.pi*A*e*1/2*density*V))**(-1)
    return W_P
def roc_constraint(eta,alpha_p,ROC,Cd0,density,A,e,W_S,beta,N_e,y):
    if y == 1: #One engine inoperative case
        W_P = ((N_e - 1)/N_e)*eta * (alpha_p/beta) *(ROC + ((4*Cd0**(1/4))/(3*np.pi*A*e)**(3/4) * np.sqrt(beta * W_S * (2/density))))**(-1)
    else:
        W_P = eta * (alpha_p/beta) *(ROC + ((4*Cd0**(1/4))/(3*np.pi*A*e)**(3/4) * np.sqrt(beta * W_S * (2/density))))**(-1)
    return W_P
def climb_gradient_constraint(eta,alpha_p,ROC_V,CD,CL,density,W_S,beta,N_e,y):
    if y == 1:
        W_P = ((N_e - 1) / N_e) * eta * (beta / alpha_p) * (1 / (ROC_V + (CD / CL))) * np.sqrt((density / 2) * ((CL) / (beta * W_S)))
    else:
        W_P = eta_prop * (beta/alpha_p) * (1/(ROC_V + (CD/CL))) * np.sqrt((density/2)*((CL)/(beta*W_S)))
    return W_P
def takeoff_constraint(alpha_p,L_to,density,h_2, k_t,N_e,y):
    if y == 1:
        W_P = alpha_p  * (1.15 * np.sqrt(((N_e) / (N_e-1))*(W_S/(L_to * k_t * density * g * np.pi * A * e))) + ((N_e) / (N_e-1))*(4*h_2/L_to)) **(-1) * np.sqrt((CL_2/W_S) *(density)/2)
    else:
        W_P = alpha_p * ((1.15 * np.sqrt((W_S / (L_to * k_t * density * g * np.pi * A * e)))) + (4 * h_2 / L_to)) ** (-1) * np.sqrt((CL_2 / W_S) * ((density) / 2))
    return W_P

propulsion_type = 2 #1. Hydrogen Combustion, 2.Hydrogen Fuel Cell, 3. Hybrid Series, 4. Hybrid Parallel Series

e_take_off = oswald_efficiency(15)
e_landing = oswald_efficiency(35)
Cd0_take_off = CD_0(15, 1)
Cd0_landing = CD_0(35, 1)
red_CD = 0.9
CD_to_props = (Cd0_take_off + (CL_max_takeoff**2 /(np.pi * A* e_take_off))) * red_CD
CD_land_props = (Cd0_landing + (CL_max_landing**2 /(np.pi * A* e_landing))) * red_CD
CD_to = Cd0_take_off + (CL_max_takeoff**2 /(np.pi * A* e_take_off))
CD_land = Cd0_landing + (CL_max_landing**2 /(np.pi * A* e_landing))


if propulsion_type == 1:
    # Approach Constraint
    W_S_approach = V_approach_constraint(rho_1524,V_approach_stall,CL_max_landing,beta_V_app_hc)
    W_S_land = s_land_constraint(s_landing_1524,C_LFL,rho_1524,CL_max_landing,beta_s_land_hc)
    W_P_cruise = cruise_contraint(eta_prop,Power_lapse(rho_cruise,rho_0), Cd0, rho_cruise, V_cruise, W_S, A, e, beta_cruise_hc)
    W_P_ROC = roc_constraint(eta_prop,Power_lapse(rho_cruise,rho_0),ROC,Cd0,rho_1524,A,e,W_S,beta_ROC_hc,2,2)
    #W_P_ROC_OEI = roc_constraint(eta_prop, Power_lapse(rho_1524, rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC_hc,2,1)
    W_P_CV = climb_gradient_constraint(eta_prop,Power_lapse(rho_1524,rho_0),ROC_V,CD_to,CL_max_takeoff,rho_1524,W_S,beta_cv,2,2)
    W_P_CV_OEI = climb_gradient_constraint(eta_prop, Power_lapse(rho_1524, rho_0), ROC_V_OEI, CD_to, CL_max_takeoff, rho_1524, W_S,beta_cv,2,1)
    W_P_TOP = takeoff_constraint(Power_lapse(rho_1524,rho_0),s_takeoff_1524,rho_1524,h2,k_t,2,2)
    W_P_TOP_OEI = takeoff_constraint(Power_lapse(rho_1524, rho_0), s_takeoff_1524, rho_1524, h2, k_t, 2, 1)
elif propulsion_type == 2:
    W_S_approach = V_approach_constraint(rho_1524, V_approach_stall, CL_max_landing, beta_V_app_fc)
    W_S_land = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_max_landing, beta_s_land_fc)
    W_P_cruise = cruise_contraint(eta_prop,alpha_p_em, Cd0, rho_cruise, V_cruise, W_S, A, e,beta_cruise_fc)
    W_P_ROC = roc_constraint(eta_prop, alpha_p_em, ROC, Cd0, rho_1524, A, e, W_S, beta_ROC_fc, 2, 2)
    #W_P_ROC_OEI = roc_constraint(eta_prop, alpha_p_em, ROC, Cd0, rho_1524, A, e, W_S, beta_ROC_fc, 2, 1)
    W_P_CV = climb_gradient_constraint(eta_prop, alpha_p_em, ROC_V, CD_to_props, CL_max_takeoff, rho_1524, W_S,beta_em, 2, 2)
    W_P_CV_OEI = climb_gradient_constraint(eta_prop, alpha_p_em, ROC_V_OEI, CD_to_props, CL_max_takeoff, rho_1524,W_S, beta_cv, 2, 1)
    W_P_TOP = takeoff_constraint(alpha_p_em, s_takeoff_1524, rho_1524, h2, k_t, 2, 2)
    W_P_TOP_OEI = takeoff_constraint(alpha_p_em, s_takeoff_1524, rho_1524, h2, k_t, 2, 1)
elif propulsion_type == 3:
    W_S_approach = V_approach_constraint(rho_1524, V_approach_stall, CL_max_landing, beta_V_app_ker)
    #W_S_approach_ice = V_approach_constraint(rho_1524, V_approach_stall, CL_max_landing, beta_V_app)
    W_S_land = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_max_landing, beta_s_land_ker)
    #W_S_land_ice = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_max_landing, beta_s_land)
    W_P_cruise = cruise_contraint(eta_prop,alpha_p_em, Cd0, rho_cruise, V_cruise, W_S, A, e, beta_cruise_ker)
    #W_P_cruise_ice = cruise_contraint(eta_prop,Power_lapse(rho_cruise,rho_0), Cd0, rho_cruise, V_cruise, W_S, A, e, beta_cruise)
    W_P_ROC = roc_constraint(eta_prop, alpha_p_em, ROC, Cd0, rho_1524, A, e, W_S, beta_cruise_ker,2,2)
    #W_P_ROC_OEI = roc_constraint(eta_prop, alpha_p_em, ROC, Cd0, rho_1524, A, e, W_S, beta_em,2,1)
    #W_P_ROC_ice = roc_constraint(eta_prop, Power_lapse(rho_1524, rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC, 2)
    #W_P_ROC_OEI_ice = roc_constraint(eta_prop, Power_lapse(rho_1524, rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC, 1)
    W_P_CV = climb_gradient_constraint(eta_prop, alpha_p_em, ROC_V, CD_to_props, CL_max_takeoff, rho_1524, W_S, beta_em, 2, 2)
    W_P_CV_OEI = climb_gradient_constraint(eta_prop, alpha_p_em, ROC_V_OEI, CD_to_props, CL_max_takeoff, rho_1524, W_S, beta_em, 2, 1)
    W_P_TOP = takeoff_constraint(alpha_p_em, s_takeoff_1524, rho_1524, h2, k_t, 2, 2)
    W_P_TOP_OEI = takeoff_constraint(alpha_p_em, s_takeoff_1524, rho_1524, h2, k_t, 2, 1)
elif propulsion_type == 4:
    W_S_approach = V_approach_constraint(rho_1524, V_approach_stall, CL_max_landing, beta_V_app_ker)
    W_S_land = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_max_landing, beta_s_land_ker)
    W_P_cruise = cruise_contraint(eta_prop, Power_lapse(rho_cruise,rho_0), Cd0, rho_cruise, V_cruise, W_S, A, e, beta_cruise_ker)
    W_P_ROC = roc_constraint(eta_prop, Power_lapse(rho_cruise,rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC_fc, 2, 2)
    W_P_CV = climb_gradient_constraint(eta_prop, Power_lapse(rho_1524,rho_0), ROC_V, CD_to_props, CL_max_takeoff, rho_1524, W_S, beta_em, 2, 2)
    W_P_CV_OEI = climb_gradient_constraint(eta_prop, Power_lapse(rho_1524,rho_0), ROC_V_OEI, CD_to_props, CL_max_takeoff, rho_1524, W_S, beta_em,2, 1)
    W_P_TOP = takeoff_constraint(Power_lapse(rho_1524,rho_0), s_takeoff_1524, rho_1524, h2, k_t, 2, 2)
    W_P_TOP_OEI = takeoff_constraint(Power_lapse(rho_1524,rho_0), s_takeoff_1524, rho_1524, h2, k_t, 2, 1)
'''
plt.vlines(W_S_approach,0,100,'b',label="Approach Speed Constraint")
plt.plot(W_S,W_P_TOP,'g',label = "Takeoff Constraint")
plt.plot(W_S,W_P_TOP_OEI,'r',label = "Takeoff Constraint (OEI)")
plt.vlines(W_S_land,0,100,'c',label ="Landing Constraint")
plt.plot(W_S,W_P_cruise,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'y',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'k',label = "Climb Gradient Constraint")
plt.plot(W_S,W_P_CV_OEI,'indigo',label = "Climb Gradient Constraint (OEI)")
plt.xlim(0,6000)
plt.ylim(0,0.5)
plt.xticks(np.arange(0,6001,500))
plt.yticks(np.arange(0,0.5,0.05))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
if propulsion_type == 1:
    plt.title("Hydrogen Combustion")
    plt.fill_between(W_S, W_P_TOP_OEI, 1, color="red", alpha=0.1)
    plt.axvspan(4282,6000,color = "red", alpha = 0.1)
    plt.fill_between(W_S, W_P_cruise, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_CV_OEI, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_ROC, 1, color="red", alpha=0.1)
    plt.plot(3532,0.0535,'o',label = "Design point")
elif propulsion_type == 2:
    plt.title("Hydrogen Fuel cell")
    plt.fill_between(W_S, W_P_TOP_OEI, 1, color="red", alpha=0.1)
    plt.axvspan(4198, 6000, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_cruise, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_CV_OEI, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_ROC, 1, color="red", alpha=0.1)
    plt.plot(3560, 0.076, 'o',label = "Design point")
elif propulsion_type == 3:
    plt.title("Hybrid Electric Series")
    plt.fill_between(W_S, W_P_TOP_OEI, 1, color="red", alpha=0.1)
    plt.axvspan(4430, 6000, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_cruise, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_CV_OEI, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_ROC, 1, color="red", alpha=0.1)
    plt.plot(3564, 0.0764, 'o',label = "Design point")
elif propulsion_type ==4:
    plt.title("Hybrid Electric Parallel Series")
    plt.fill_between(W_S, W_P_TOP_OEI, 1, color="red", alpha=0.1)
    plt.axvspan(4430, 6000, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_cruise, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_CV_OEI, 1, color="red", alpha=0.1)
    plt.fill_between(W_S, W_P_ROC, 1, color="red", alpha=0.1)
    plt.plot(3554, 0.07, 'o',label = "Design point")
plt.legend(loc = "upper right")
plt.grid()
plt.show()

print("Hydrogen Combustion W/S = 3552")
print("Hydrogen Combustion W/P =  0.0535")
print("Fuel Cell W/S = 3560")
print("Fuel Cell W/P = 0.0763 ")
print("Hybrid Series W/S = 3564")
print("Hybrid Series W/P =  0.0764")
print("Hybrid Parallel Series W/S = 3554")
print("Hybrid Parallel Series W/P = 0.07")'''
W_P_design = 0.0763
W_S_design = 3560