import numpy as np
from matplotlib import pyplot as plt
# Take off Distance Constraint
g= 9.80665                      #m/s^2 Gravitational
lambda_trop = -6.5/1000         #
R = 287.0528                    #
rho_0 = 1.225                   #kg/m^3
s_landing = 1370                #m
eff_prop = 0.85                 #- (turboprop)
V_cruise = 275 * 0.51444444     #m/s
h = 280*100*0.3048              #m
W_S = np.arange(1,5000,1)       #N/m^2

#CONSTRAINTS DEPENDING ON THE DESIGN
A = 12                          #- (12-14)
e = 0.7                         #-
CL_max = 1.8                    #- (1.2-1.8 twin engine) & (1.5-1.9 turboprop) max lift coefficient
Cd0 = 0.0335                    #- calculated in Class I conventional: Cdo =Swet_S*Cfo
ROC = 20                        #m/s higher more constraints
ROC_V = 0.083                   #- gradient
CL_to_max = 2.5                 #- (1.9-3.3)
CL_land = 1.9                   #- ? (1.7-2.1) (max?)
TOP = 300                       #- (420-460) -> from Raymer graph
V_stall = 61 * 0.514444         #m/s maximum 61kts for CS23

#INTERMEDIATE CALCULATIONS
CL_to = CL_to_max/1.1**2
V_to = 1.1 * V_stall
CD_to = Cd0 + (CL_to**2 /(np.pi * A* e))
rho1524_rho0 = (1 +((lambda_trop* 1524)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_rho0 = (1 +((lambda_trop* h)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_rho0 * rho_0
W_S_stall = 1/2 * rho_0 * V_stall**2 * CL_max

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to * rho_rho0 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_0 * s_landing/0.5915)/(2*0.95) #Change to CS25 regulation
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_0))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + CD_to/CL_to)*(np.sqrt((2/rho_0)*(1/CL_to))))

#plt.vlines(W_S_stall,0,0.4,'b',label="V_stall")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,0.4,'k',label ="Landing")
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,1500)
plt.ylim(0,0.4)
plt.xticks(np.arange(0,5000,500))
plt.yticks(np.arange(0,0.41,0.05))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.legend(loc = "upper right")
plt.grid()
plt.show()