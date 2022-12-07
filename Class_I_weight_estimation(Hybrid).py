import numpy as np
from matplotlib import pyplot as plt
#Constants
g = 9.80665
lambda_trop = -6.5/1000
R = 287.0528
R_norm = 1000 * 1852            #Range in meters
E = 45 * 60                     #Loiter endurance in seconds
V_cruise = 275 * 0.51444444     #knt -> m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
R_div = 200*1000                #m (TBD)  --> to be determined in literature
f_con = 5/100                   #-
s_takeoff_1524 = 1370           #Takeoff Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
s_landing_1524 = 1370           #Landing Distance at 1524 m above mean sea-level (ISA + 10 degree) (m)
rho_0= 1.225                    #ISA + 10 ◦C day (kg/m3) ADD TEMPERATURE DEVIATION
rho_1524= 1.01893               #1524m ISA + 10 ◦C day (kg/m3)
rho_1524_rho0 = rho_1524/rho_0
rho_cruise_rho0 = (1 +((lambda_trop* h_cruise)/(288.150))) ** (-1*(g/(R*lambda_trop)+1))
rho_cruise = rho_cruise_rho0 * rho_0
eff_prop = 0.85              #Change with Literature
m_payload = 5728
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
TOP = 430                          #Change with Literature Reference to slide (420-460) -> from Raymer graph
ROC = 10.2                         #Change with CS25 and literature or Requirement (Rate of Climb)
ROC_V = 0.0024                     #Change with CS25 and literature or Requirement (Climb Gradient) ROC/V
V_stall = 54                       #Change with CS25 or Requirement
a = 0.5464
b = 1439*g

#CALCULATIONS for the GRAPHS
W_P_TOP = TOP/ (W_S) * CL_to * rho_1524_rho0 #(Use this if it's at a different altitude then sea level)
# Landing Distance Constraint
W_S_land = (CL_land * rho_1524 * s_landing_1524/0.5847)/(2*0.95)
# Cruise Speed Constraint
W_P_cru = eff_prop * (rho_cruise_rho0)**(3/4) * ((((Cd0*1/2*rho_cruise*V_cruise**3)/W_S)+(W_S/(np.pi*A*e*1/2*rho_cruise*V_cruise)))**(-1))
# Rate of Climb Constraint
W_P_ROC = eff_prop / (ROC + ((np.sqrt(W_S)*np.sqrt(2/rho_1524))/(1.345*((A*e)**(3/4))/(Cd0**(1/4)))))
# Climb Gradent Constraint
W_P_CV = eff_prop / (np.sqrt(W_S)*(ROC_V + (CD_to/CL_to))*(np.sqrt((2/rho_1524)*(1/CL_to))))
#Stall Constraint
W_S_stall = 1/2 * rho_1524 * V_stall**2 * CL_to

'''plt.vlines(W_S_stall,0,100,'b',label="V_stall")
plt.plot(W_S,W_P_TOP,'r',label = "Takeoff Constraint")
plt.vlines(W_S_land,0,100,'k',label ="Landing")
plt.plot(W_S,W_P_cru,'m',label = "Cruise Constraint")
plt.plot(W_S,W_P_ROC,'c',label = "Rate of Climb Constraint")
plt.plot(W_S,W_P_CV,'y',label = "Climb Gradient Constraint")
plt.xlim(0,5000)
plt.ylim(0,0.5)
plt.xticks(np.arange(0,5001,500))
plt.yticks(np.arange(0,0.5,0.05))
plt.xlabel("W/S (N/m^2)")
plt.ylabel("W/P (N/W)")
plt.legend(loc = "upper right")
plt.grid()
plt.show()'''

#Mass Preliminary Calculation
W_P_design = 0.0468
W_S_design = 3112
m_turboprop = 1074.5
m_em = 50
MTOW_design = 207146                    #N
P_max = MTOW_design / W_P_design
S = MTOW_design / W_S_design
#print(P_max/1000, "kW Max Power")
#print(S,"m^2 Surface Area ")
#Degree of Hybridization of Energy (He) *Could be defined by each split point or total journey
t_cruise_500 = 4476
t_cruise_full = 17568
t_climb1 = 222
t_climb2 = 600
t_climb3 = 975
t_descent1 = 503
t_descent2 = 377
L_D_cruise = 16.7
L_D_to = CL_to/CD_to
L_D_land = 6.1
t_total_500 = t_cruise_500 + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
t_total_full = t_cruise_full + t_climb1 + t_climb2 + t_climb3 + t_descent1 + t_descent2
delta_v = V_cruise
E_total_climb1 = (MTOW_design*71.7)/ (L_D_to) * t_climb1 + (MTOW_design/g * 72**2) + MTOW_design*6.8*t_climb1
E_total_climb2 = (MTOW_design*107.9)/ (L_D_cruise) * t_climb2 + (MTOW_design/g * 36**2) + MTOW_design*5.08*t_climb2
E_total_climb3 = (MTOW_design*107.9)/ (L_D_cruise) * t_climb3 + MTOW_design*4.06*t_climb3
E_total_cruise_500 = (MTOW_design*V_cruise)/ (L_D_cruise) * t_cruise_500 + (MTOW_design/g * 33.4**2)
E_total_cruise_full = (MTOW_design*V_cruise)/ (L_D_cruise) * t_cruise_full + (MTOW_design/g * 33.4**2)
E_total_descent1 = (MTOW_design*138.9)/ (L_D_cruise) * t_descent1 + (MTOW_design/g * 2.57**2) - MTOW_design*10.89*t_descent1
E_total_descent2 = (MTOW_design*103)/ (L_D_land) * t_descent2 + (MTOW_design/g * 36**2) - MTOW_design*8.08*t_descent2
E_total_500 = E_total_climb1 + E_total_climb2 + E_total_climb3 + E_total_descent1 +E_total_descent2 +E_total_cruise_500
E_total_full = E_total_climb1 + E_total_climb2 + E_total_climb3 + E_total_descent1 +E_total_descent2 +E_total_cruise_full

print(np.round(E_total_500/10**6,5),"MJ, 500Nmi")
print(np.round(E_total_full/10**6,5),"MJ, Equivalent Range")
'''print(E_total_climb1/10**6, "MJ Climb1")
print(E_total_climb2/10**6, "MJ Climb2")
print(E_total_climb3/10**6, "MJ Climb3")
print(E_total_cruise_full/10**6, "MJ Cruise 500nmi")
print(E_total_cruise_full/10**6, "MJ Cruise full")'''

E_nc = 0 * E_total_full
He = E_nc / E_total_full        #Energy of non consumable(battery) / Total Energy
#print(E_nc/10**6,"MJ for non consumable")
#print((E_total_full - E_nc)/10**6,"MJ for consumable")
eta_stt = 0.85 * 0.45       #Efficiency chain from shaft-to-thrust
eta_fuel_cell = 0.7 * 0.95 * 0.99
eta_btt = 0.95 * 0.75       #Efficiency chain from battery-to-thrust
NoD_ice = 2                 #Number of turboprop engines
NoD_em = 4                 #Number of electric motor engines
P_ice = (E_total_full - E_nc)/ (eta_stt * t_total_full * NoD_ice)
P_em = E_nc/ (eta_btt * t_total_full * NoD_em)
P_ice_climb1 = (E_total_climb1)/ (eta_fuel_cell * t_climb1 * NoD_ice)
P_ice_climb2 = (E_total_climb2)/ (eta_fuel_cell * t_climb2 * NoD_ice)
P_ice_climb3 = (E_total_climb3)/ (eta_fuel_cell * t_climb3 * NoD_ice)
P_ice_cruise_500 = (E_total_cruise_500)/ (eta_fuel_cell * t_cruise_500 * NoD_ice)
P_ice_cruise_full = (E_total_cruise_full)/ (eta_fuel_cell * t_cruise_full * NoD_ice)
print(np.round(P_ice_climb1/10**3), "kW Climb1")
print(np.round(P_ice_climb2/10**3), "kW Climb2")
print(np.round(P_ice_climb3/10**3), "kW Climb3")
print(np.round(P_ice_cruise_500/10**3), "kW Cruise 500nmi")
print(np.round(P_ice_cruise_full/10**3), "kW Cruise full")

#print(P_ice/1000,"kW for ICE per engine")
#print(P_em/1000,"kW for Electric Battery per engine")
#Degree of Hybridization of Power(Hp)
# Choice between Parallel and Series needs to be made
#If Parallel:
H_p_para = P_em*NoD_em / P_max
#print(np.round(H_p_para,2)*100,"% Hybridlization")
'''#If Series:
H_p_ser = P_em_max / P_ice_max'''
tf =  0                      #Trap fuel time step
BSFC= 0.48*(0.45/(745*3600)) #Brake-specific fuel consumption
ddp = 0.8                   #Deep discharge protection
E_bat = 2.7*10**6           #Total Battery Energy per piece
m_fuel_ice = (1+tf)*P_ice*NoD_ice*BSFC*t_total_500
m_bat = (1+ddp) * (E_nc/(eta_btt*E_bat))
m_OE = (a * MTOW_design + b)/g - m_turboprop        #Maximum Takeoff Mass
m_propulsion = m_turboprop + m_em*NoD_em
m_MTOW = m_OE + m_fuel_ice + m_payload + m_bat + m_propulsion
'''print(m_fuel_ice, "Fuel Mass(kg)")
print(m_bat,"Battery Mass(kg)")
print(m_OE, "Operational Empty mass without Engine(kg)")
print(m_propulsion, "Propulsion Mass (kg)")
print(m_MTOW, "Calculated MTOW (kg)")
print(MTOW_design/g, "MTOW Design(kg)")'''