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
g = 9.80665
R_norm = 1000 * 1852            #Range in meters
E = 45 * 60                     #Loiter endurance in seconds
V_cruise = 275 * 0.51444444     #m/s (TAS)
h_cruise = 280*100 * 0.3048     #m
R_div = 0                       #m (TBD)  --> to be determined in literature
f_con = 5/100                   #-
e_kero = 42.9                   #MJ/kg Specific Energy Kerosene
e_atj = 43.2                    #MJ/kg Specific Energy SAF(ATJ)
e_lh2 = 142                     #MJ/kg Specific Energy Liquid Hydrogen
e_bat = 1.1                     #MJ/kg Specific Energy Battery (assuming 300Wh/kg)

#Constant on estimation
PAX = 50                        #Number of passengers
WCargo =                        #Cargo Weight (TBD)
eta_EM = 0.95                   # Electric motor efficiency
PSFC =                          # Power Specific Fuel Consumption
TSFC =                          # Thrust Specific Fuel Consumption
eta_p_home = 0.7                # Propulsive efficiency - homebuilts
eta_p_single = 0.80             # Propulsive efficiency - single engine/motor props
eta_p_twin = 0.82               # Propulsive efficiency - twin engine props
eta_p_turbo = 0.85              # Propulsive efficiency - regional turboprops
a_p = 0.5464                    # turboprop     WE = a*MTOW+b
a_j = 0.4985                    # turbojet
b_p = 1439                      # turboprop
b_j = 1782.3                    # turbojet

# Variable Constants depending on aircraft features
A =                             # Aspect ratio (CHANGE) - high A for slender wing
#Cdo calculations
Psi = 0.0075 #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97   #span efficiency factor (value based on Roelof reader p.46)
e = 1/((np.pi)*A*Psi+(1/phi))
Cfe =                           #Equivalent skin friction coefficient - depending on aircraft from empirical estimation
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe * Swet_S
CL = np.sqrt(np.pi()*Cd0*A*e)
CD = 2 * Cd0

# Depending on energy source
e_f =                           # Specfic Energy per fuel type
eta_eng =                       # Engine efficiency
eta_p =                         # Propulsive efficiency depending on airplane types
a =                             #Linear regression for OEW
b =                             #Linear regression for OEW

#calculation Turboprop/piston
eta_eng = (1/e_f)*(1/PSFC)      # Thermodynamic efficiency of engine

#calculation Jet engine
eta_eng = (V_cruise/TSFC)*(1/e_f)*(1/eta_p)

R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))
Req = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + E*V_cruise

#take-off, climb, descent, deceleration
R = R_lost*(1+f_con)
mfuel_MTO = 1- np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))
mbat_MTO = R/(eta_EM*eta_p*(e_bat/g)*(CL/CD))

#cruise + loiter cruise
R = (R_norm)*(1+f_con) + 1.2 * R_div + E*V_cruise
mfuel_MTO = 1- np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))
mbat_MTO = R/(eta_EM*eta_p*(e_bat/g)*(CL/CD))

#cruise + loiter cruise + take-off, clmib, descent, deceleration
R = Req
mfuel_MTO = 1- np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))
mbat_MTO = R/(eta_EM*eta_p*(e_bat/g)*(CL/CD))

##dependent on aircraft fuel use
Mff =                                   #dependent on energy source
Mres = 0.25

WPAX = 200*0.453592*PAX
WPAXBAGGAGE = 40*0.453592*PAX


WPLtot = WPAX + WPAXBAGGAGE + WCargo
MTOW = (b + WPLtot)/(Mff-a)
WOE = a*MTOW + b
WF = (Mres+1)*(1-Mff)*MTOW                  ##TO BE CHECKED

##REWRITTEN main FORMULA & add output formulas
