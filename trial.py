import numpy as np
from matplotlib import pyplot as plt
#Constants
'''g = 9.80665
V_cruise = 55
V_stall = 25
h_cruise = 3000     #m
s_takeoff = 200           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
W_S = np.arange(1,1000,1)
A = 8
e = 0.8
Cd0 = 0.027
CL_0 = 0.25
CL_to = 2.0
ROC = 4.5
eff_prop = 0.85
k = 1/ (np.pi * A * e)
CD_to = Cd0 + k *(CL_to - CL_0)**2
rho = 1.225
rho_cruise = 0.909
friction_g = 0.72
# Take off Distance Constraint
V_to = 1.1 * np.sqrt(2 * W_S/(rho * CL_to))
T_W_TOD = V_to**2 / (2*g*s_takeoff) + ((1/2 *rho *V_to*2)*CD_to)/W_S + friction_g *(1 - (1/2 *rho *V_to*2*CL_to)/W_S)
P_W_TOD = T_W_TOD * V_to / eff_prop

# Cruise Speed Constraint
T_W_crui = ((1/2*rho_cruise*V_cruise**2)/W_S) *(Cd0 + k * ((W_S/(1/2*rho*V_cruise**2))-CL_0)**2)
P_W_crui_0 = T_W_crui * V_cruise / eff_prop
P_W_crui = P_W_crui_0 * (rho_cruise/rho - ((1-(rho_cruise/rho))/7.55))
#ROC Constraint
V_ldmax = ((W_S **2 * 1/(rho/2)**2) * (k/(Cd0 + k * CL_0**2)))**(1/4)
H = W_S * 1/(CL_to/CD_to*rho) * (2*k*Cd0*CL_to/CD_to+1)/(Cd0+k*CL_0**2)
V_y = np.sqrt(H- np.sqrt(H**2 - V_ldmax**4))
T_W_ROC = ROC/V_y+ (1/2 * rho * V_to**2)/W_S * (Cd0 + k * (W_S/(1/2 * rho * V_to**2)-CL_0)**2)
P_W_ROC = T_W_ROC * V_to / eff_prop
# Stall Constraint
W_S_stall = rho/2 * V_stall**2 * CL_to


plt.vlines(W_S_stall,0,150,'b',label = "V_stall")
plt.plot(W_S,P_W_TOD,'r',label = "Takeoff Constraint")
plt.plot(W_S,P_W_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,P_W_crui,'m',label = "Cruise Constraint")
plt.xlabel("W/S")
plt.ylabel("P/W")
plt.xlim(0,900)
plt.ylim(0,150)
plt.legend(loc = "upper right")
plt.show()'''
# Take off Distance Constraint
g= 9.80665
lambda_trop = -6.5/1000
R = 287.0528
V_stall = 61 * 0.514444
CL_max = 1.8
rho_0 = 1.225
V_to = 1.1 * V_stall
CL_to = 1.3
CL_land = 1.6
TOP = 58.18
s_landing = 550
eff_prop = 0.8
A = 6
e = 0.70
Cd0 = 0.0335
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
V_cruise = 70
h = 2000
rho_rho0 = (1 +((lambda_trop* h)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_2000 = rho_rho0 * rho_0
ROC = 5
ROC_V = 0.083
W_S = np.arange(1,1500,1)
W_S_stall = 1/2 * rho_0 * V_stall**2 * CL_max

W_P_TOP = TOP/ (W_S) * CL_to #* rho_rho0 (Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_0 * s_landing/0.5915)/(2*0.95) #Change to CS25 regulation
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_rho0)**(3/4) * ((((Cd0*1/2*rho_2000*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_2000*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + CD_to/CL_to)*(np.sqrt((2/rho_0)*(1/CL_to))))
plt.vlines(W_S_stall,0,0.4,'b',label="V_stall")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,0.4,'k',label ="Landing")
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,1500)
plt.ylim(0,0.4)
plt.xticks(np.arange(0,1501,500))
plt.yticks(np.arange(0,0.41,0.05))
plt.legend(loc = "upper right")
plt.grid()
plt.show()