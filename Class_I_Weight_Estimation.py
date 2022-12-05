import numpy as np

# Outputs: MTOW OEW and WF
# Relationships: MTOW = WOE + WF + WPL
# Relationships: MTOW = WE + WF + WPLtot + Wtfo
# LOOK INTO MF/MTO FORMULA FROM ROLOEF READER
# LIQUID FUELS WILL HAVE THE SAME FORMULA AS JET-A FUEL
# BATTERIES WILL HAVE FORMULA WITHOUT ln()

#Handy unit conversion
#1 lbs = 0.453592 kg
#1 ft = 0.3048 m

#Constants
g = 9.80665                     #gravity
R_norm = 1852000                #Range in m
E = 45 * 60                     #Loiter endurance in seconds
V_cruise = 275 * 0.51444444     #m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
R_div = 200000                  #m
f_con = 5/100                   #-
PAX = 50                        #Number of passengers
Mres = 0.25                     #-

#Constants depending on the aircraft
e_kero = 42.9                   #MJ/kg Specific Energy Kerosene
e_atj = 43.2                    #MJ/kg Specific Energy SAF(ATJ)
e_lh2 = 142                     #MJ/kg Specific Energy Liquid Hydrogen
e_bat = 1.1                     #MJ/kg Specific Energy Battery (assuming 300Wh/kg)
PSFC = 0.48*(0.45/(745*3600))   #kg/N/s Power Specific Fuel Consumption    ---> to be checked
TSFC = 0.67*(0.45/(745*3600))   #kg/N/s Thrust Specific Fuel Consumption   ---> to be checked
eta_EM = 0.95                   #- electric motor efficiency
eta_p_single = 0.80             #- propulsive efficiency - single engine/motor props
eta_p_twin = 0.82               #- propulsive efficiency - twin engine props
eta_p_turbo = 0.85              #- propulsive efficiency - regional turboprops
a_p = 0.5464                    #- turboprop     WE = a*MTOW+b
a_j = 0.4985                    #- turbojet
b_p = 1439                      #N turboprop
b_j = 1782.3                    #N turbojet
Psi = 0.0075                    #- parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #- span efficiency factor (value based on Roelof reader p.46)

#TO BE CHANGED DEPENDING ON THE DESIGN
eta_p =eta_p_turbo              #- propulsive efficiency depending on airplane types (TYPES: single, twin, regional turboprop)
e_f = e_lh2*1000000             #J/kg specfic Energy per fuel type  (TYPES: kerosene, SAF, LH)

#Intermediate calculations
eta_eng_prop = (1/e_f)*(1/PSFC)                     # thermodynamic efficiency of engine turbo/piston
eta_eng_jet = (V_cruise/TSFC)*(1/e_f)*(1/eta_p)     # thermodynamic efficiency of engine jet

#TO BE CHANGED DEPENDING ON THE DESIGN
A = 12                          #- aspect ratio -> high A for slender wing
Cfe = 0.003                     #- equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6                      #- (6.0-6.2) wetted area ratios -> depending on airframe structure
eta_eng = eta_eng_prop          #- engine efficiency             (TYPES: jet, propeller)
a =    a_p                      #- linear regression for OEW     (TYPES: turboprop, turbojet)
b =       b_p                   #N linear regression for OEW     (TYPES: turboprop, turbojet)

#Calculations Aerodynamic Characteristics
e = 1/((np.pi)*A*Psi+(1/phi))   #-
Cd0 = Cfe * Swet_S              #-
CL = np.sqrt(np.pi*Cd0*A*e)     #-
CD = 2 * Cd0                    #-

#Calculations Range
R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))                #m
R = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + (E*V_cruise)                  #m

#Calculations fuel fraction
mfuel_MTO_FULL = 1-np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))                 #- fuel
mbat_MTO_FULL = R/(eta_EM*eta_p*(e_bat/g)*(CL/CD))                            #- battery

#TO BE CHANGED DEPENDING ON THE DESIGN
Mused = mfuel_MTO_FULL                                                        #-

WPAX = 200*0.453592*PAX*g + 3*190*0.453592*g                                  #N
WPAXBAGGAGE = 40*0.453592*PAX*g +2*30*0.453592*g                              #N Crew is bagageless


#OUTPUTS
WPLtot = WPAX + WPAXBAGGAGE                                                   #N
MTOW = (b + WPLtot)/(1-a-Mused*(1+Mres))                                      #N
WOE = a*MTOW + b                                                              #N
WF = MTOW*(Mused*(1+Mres))                                                    #N



print("MTOW =", MTOW)
print("Fuel weight =", WF)
print("Operational Empty Weight =", WOE)
#print("R =",R)
#print("eta_engine = ", eta_eng)
#print("eta_p", eta_p)
#print("e_f = ",e_f)
print("CL =",CL)
print("CD =", CD)
#print("Rlost =", R_lost, "[m]")
#print("PSFC =", PSFC)
#print("Wpl =", WPLtot, "[N]")
