import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
lambda_trop = -6.5/1000
R = 287.0528
kts = 0.51444444
ft = 0.3048
lbs = 0.453592
V_cruise =  70 #275 * kts     #kts -> m/s (TAS)
h_cruise = 1800  #280*100 * ft     #m
s_takeoff_1524 = 750#4500 * ft           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 750 #1/ 0.6 * (4500 * ft)           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
rho_0= 1.225                    #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524= 1.01893               #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
eff_prop = 0.80 #0.85                 #Change with Literature
PAX = 50
WPAX = 200*lbs*PAX*g                                 #N
WPAXBAGGAGE = 40*lbs*PAX*g                           #N Crew is bagageless
m_payload = (WPAX + WPAXBAGGAGE) / g
W_S = np.arange(1,4000,1)
print(rho_cruise)

##Cdo calculations
Psi = 0.0075                    #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #span efficiency factor (value based on Roelof reader p.46)
A = 9 #12                           #Aspect Ratio (12-14) #Reference to ATR 72
e = 0.75 #1/(np.pi*A*Psi+(1/phi))
Cfe = 0.0045                     #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1                     #(6.0-6.2) wetted area ratios -> depending on airframe structure
Cd0 = 0.049 #Cfe * Swet_S
#Aerodynamic Estimations
CL_max = 1.4 #1.9                       #(1.2-1.8 twin engine) & (1.5-1.9 turboprop) max lift coefficient
CL_to = 1.7#2.1                        #Change with Estimate (1.7-2.1)
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
CL_land = 2.1 #2.6                      #Change with Estimate (1.9-3.3)
CD_land = Cd0 + (CL_land**2 /(np.pi * A* e))
TOP = 430                          #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC = 4                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0032
ROC_V_OEI = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_approach = 141 * 0.514444         #Requirement
V_approach_stall = V_approach /1.23  #CS 25 requirement of V_stall_land = V_approach / 1.23
a = 0.5088
b = 1199.7
beta_V_app = 0.85
beta_s_land = 0.85
beta_cruise = 0.95
beta_ROC = 0.95
beta_cv = 1
beta_s_to = 1
beta_em = 1
C_LFL = 0.6 #0.45                        #Landing field length coefficient s^2/m
alpha_p_em = 1
CL_2 = 1.18 #(1.13)**2 * CL_to
k_t = 0.85
h2 = 50 * ft
V_s0 = 31
def Power_lapse(rho,rho_0):
    alpha_p_ice = (rho/rho_0) ** (3/4)
    return alpha_p_ice
def V_approach_constraint(density,V_approach_stall,CL_land,beta):
    W_S = 1/2 * density * V_approach_stall **2 * CL_land * 1/beta
    return W_S
def s_land_constraint(s_land,C_LFL,density,CL_land,beta):
    W_S = (s_land / C_LFL) * (density * CL_land / 2) * 1/beta
    return W_S
def cruise_contraint(eff_prop,alpha_p,Cd0,density,V,W_S,A,e,beta):
    W_P = eff_prop * (alpha_p/beta) *((Cd0*1/2*density*V**3)/(beta*W_S)+ (beta*W_S)/(np.pi*A*e*1/2*density*V))**(-1)
    return W_P
def roc_constraint(eff_prop,alpha_p,ROC,Cd0,density,A,e,W_S,beta,N_e,y):
    if y == 1: #One engine inoperative case
        W_P = ((N_e - 1)/N_e)*eff_prop * (alpha_p/beta) *(ROC + ((4*Cd0**(1/4))/(3*np.pi*A*e)**(3/4) * np.sqrt(beta * W_S * (2/density))))**(-1)
    else:
        W_P = eff_prop * (alpha_p/beta) *(ROC + ((4*Cd0**(1/4))/(3*np.pi*A*e)**(3/4) * np.sqrt(beta * W_S * (2/density))))**(-1)
    return W_P
def climb_gradient_constraint(eff_prop,alpha_p,ROC_V,CD,CL,density,W_S,beta,N_e,y):
    if y == 1:
        W_P = ((N_e - 1) / N_e) * eff_prop * (beta / alpha_p) * (1 / (ROC_V + (CD / CL))) * np.sqrt((density / 2) * ((CL) / (beta * W_S)))
    else:
        W_P = eff_prop * (beta/alpha_p) * (1/(ROC_V + (CD/CL))) * np.sqrt((density/2)*((CL)/(beta*W_S)))
    return W_P
def takeoff_constraint(alpha_p,L_to,density,h_2, k_t,N_e,y):
    if y == 1:
        W_P = alpha_p  * (1.15 * np.sqrt(((N_e - 1) / N_e)*(W_S/(L_to * k_t * density * g * np.pi * A * e))) + ((N_e - 1) / N_e)*(4*h_2/L_to)) **(-1) * np.sqrt((CL_2/W_S) *(density)/2)
    else:
        W_P = alpha_p * (1.15 * np.sqrt((W_S / (L_to * k_t * density * g * np.pi * A * e))) + (4 * h_2 / L_to)) ** (-1) * np.sqrt((CL_2 / W_S) * (density) / 2)
    return W_P

propulsion_type = 1 #1. Hydrogen Combustion, 2.Hydrogen Fuel Cell, 3. Hybrid Series, 4. Hybrid Parallel Series

'''
if propulsion_type == 1 or 2:
    # Approach Constraint
    W_S_approach = V_approach_constraint(rho_1524,V_approach_stall,CL_land,beta_V_app)
    W_S_land = s_land_constraint(s_landing_1524,C_LFL,rho_1524,CL_land,beta_s_land)
    W_P_cruise = Cruise_contraint(eff_prop,Power_lapse(rho_cruise,rho_0), Cd0, rho_cruise, V_cruise, W_S, A, e, beta_cruise)
    W_P_ROC = ROC_constraint(eff_prop,Power_lapse(rho_1524,rho_0),ROC,Cd0,rho_1524,A,e,W_S,beta_ROC,2)
    W_P_ROC_OEI = ROC_constraint(eff_prop, Power_lapse(rho_1524, rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC, 1)
    W_P_CV = climb_gradient_constraint(eff_prop,Power_lapse(rho_1524,rho_0),ROC_V,CD_to,CL_to,rho_1524,W_S,beta_cv,2)
    W_P_CV_OEI = climb_gradient_constraint(eff_prop, Power_lapse(rho_1524, rho_0), ROC_V_OEI, CD_to, CL_to, rho_1524, W_S,beta_cv, 1)
if propulsion_type == 3 or 4:
    # Approach Constraint
    W_S_approach_em = V_approach_constraint(rho_1524, V_approach_stall, CL_land, beta_em)
    W_S_approach_ice = V_approach_constraint(rho_1524, V_approach_stall, CL_land, beta_V_app)
    W_S_land_em = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_land, beta_em)
    W_S_land_ice = s_land_constraint(s_landing_1524, C_LFL, rho_1524, CL_land, beta_s_land)
    W_P_cruise_em = Cruise_contraint(eff_prop,alpha_p_em, Cd0, rho_cruise, V_cruise, W_S, A, e, beta_em)
    W_P_cruise_ice = Cruise_contraint(eff_prop,Power_lapse(rho_cruise,rho_0), Cd0, rho_cruise, V_cruise, W_S, A, e, beta_cruise)
    W_P_ROC_em = ROC_constraint(eff_prop, alpha_p_em, ROC, Cd0, rho_1524, A, e, W_S, beta_em, 2)
    W_P_ROC_OEI_em = ROC_constraint(eff_prop, alpha_p_em, ROC, Cd0, rho_1524, A, e, W_S, beta_em, 1)
    W_P_ROC_ice = ROC_constraint(eff_prop, Power_lapse(rho_1524, rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC, 2)
    W_P_ROC_OEI_ice = ROC_constraint(eff_prop, Power_lapse(rho_1524, rho_0), ROC, Cd0, rho_1524, A, e, W_S, beta_ROC, 1)
'''
W_S_approach = rho_0/ 2 * V_s0 **2 * CL_land
W_P_TOP = takeoff_constraint(1,s_takeoff_1524,rho_0,h2,k_t,1,2)
W_S_land = s_land_constraint(s_landing_1524,C_LFL,rho_0,CL_land,beta_em)
W_P_cru = cruise_contraint(eff_prop,1,0.026,1.009,70,W_S,9,0.71,beta_em)
W_P_ROC = roc_constraint(0.8,1,2.0,0.026,1.23,9.0,0.71,W_S,1,1,2)
W_P_CV = climb_gradient_constraint(0.8,1,0.083,0.16,1.5,1.23,W_S,1,1,2)
W_P_TOP = takeoff_constraint(1,750,1.23,50*ft,k_t,1,2)
plt.vlines(W_S_approach,0,100,'b',label="Approach Speed Constraint")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,100,'k',label ="Landing Constraint")
#plt.axvspan(3171,4000,color = "red", alpha = 0.1)
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
#plt.fill_between(W_S,W_P_cru,1,color = "red",alpha = 0.1)
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,2000)
plt.ylim(0,0.5)
plt.xticks(np.arange(0,2001,500))
plt.yticks(np.arange(0,0.51,0.1))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.legend(loc = "upper right")
plt.grid()
plt.show()

'''#Mass Preliminary Calculation
W_P_design = 0.0467
W_S_design = 3171
MTOW_design = 20281 * g                  #N
#print(MTOW_design/W_S_design,"m^2")
#print(MTOW_design/W_P_design/10**3,"kW")'''
print(W_P_CV[501])
print(W_P_CV[1001])
print(W_P_CV[1501])
print(W_P_CV[2001])