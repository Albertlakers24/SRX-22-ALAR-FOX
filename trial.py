import numpy as np
from matplotlib import pyplot as plt
# Take off Distance Constraint
g= 9.80665
lambda_trop = -6.5/1000
R = 287.0528
V_stall = 61 * 0.514444     #Check with CS25
CL_max = 1.8                #Change with Estimate
rho_0 = 1.225
V_to = 1.1 * V_stall
CL_to = 1.3                 #Change with Estimate
CL_land = 1.6               #Change with Estimate
TOP = 58.18                 #Change with Literature
s_landing = 550             #Change with Requirement
eff_prop = 0.8              #Change with Literature
e = 0.70
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
V_cruise = 70
h = 2000
rho_rho0 = (1 +((lambda_trop* h)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_2000 = rho_rho0 * rho_0
ROC = 5
ROC_V = 0.083
W_S = np.arange(1,1500,1)
#Stall Constraint
W_S_stall = 1/2 * rho_0 * V_stall**2 * CL_max
#Takeoff Constraint
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