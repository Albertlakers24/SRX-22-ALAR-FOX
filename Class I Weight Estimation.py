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
R_div = 0                       #m (CHANGE)
f_con = 5/100                   #-
e_kero = 42.9                   #MJ/kg Specific Energy Kerosene
e_atj = 43.2                    #MJ/kg Specific Energy SAF(ATJ)
e_lh2 = 142                     #MJ/kg Specific Energy Liquid Hydrogen
e_bat = 1.1                     #MJ/kg Specific Energy Battery (assuming 300Wh/kg)


#Constant on estimation
A =                             #Aspect ratio (CHANGE) - high A for slender wing
e =                             #Oswald factor (CHANGE)

PAX = 50        #Number of passengers
#WE = ##linear regression relating to MTOW depending on the aircraft
OEW = ## linear regression relating to MTOW dependiin on the aircraft
a = 0.5422 or 0.4985            #Linear regression from MTOW/ OEW (turboprop or turbrojet)
b = 1455.2 or 1782.3            #Linear regression from MTOW/ OEW (turboprop or turbrojet)
eta_p =                         #propeller efficiency   -> maximize (CHANGE)
c_p =                           #propeller              -> minimize (CHANGE)
c_j =                           #jet                    -> minimize (CHANGE)
se =                            #Specific Energy (MJ/kg) (CHANGE)
eta_EM =                        # Electric motor efficiency
eta_p =                         # Propulsive efficiency
eta_eng =                       # Thermodynamic efficiency of engine

#linear regression relating to MTOW depending on the aircraft   WE = a*MTOW+b
a = 0.5464 or 0.4985            #Linear regression from MTOW/ OEW (turboprop or turbojet)
b = 1439 or 1782.3            #Linear regression from MTOW/ OEW (turboprop or turbojet)

#Wtfo = Mtfo * MTOW    (trapped fuel and oil)
#WF    use the fuel fraction method
WF = Mused*MTOW*(Mres + 1)

#Cdo calculations
Psi = 0.0075 #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97   #span efficiency factor (value based on Roelof reader p.46)
e = 1/((np.pi)*A*Psi+(1/phi))
Cfe =                           #Equivalent skin friction coefficient - depending on aircraft from empirical estimation
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe*Swet_S

#W4/W5
CL = np.sqrt(np.pi()*Cd0*A*e)
CD = 2 * Cd0
R_lost = 1 / 0.7 * (CL/CD) *(h_cruise + (V_cruise**2 / (2*g)))
Req = (R_norm + R_lost)*(1+f_con) + 1.2 * R_div + E*V_cruise

#take-off, climb, descent, deceleration
R = R_lost*(1+f_con)
mfuel_MTO = 1- np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))
mbat_MTO = R/(eta_EM*eta_p*(e_f/g)*(CL/CD))

#cruise + loiter cruise
R = (R_norm)*(1+f_con) + 1.2 * R_div + E*V_cruise
mfuel_MTO = 1- np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))
mbat_MTO = R/(eta_EM*eta_p*(e_f/g)*(CL/CD))

#cruise + loiter cruise + take-off, clmib, descent, deceleration
R = Req
mfuel_MTO = 1- np.exp(-R/(eta_eng*eta_p*(e_f/g)*(CL/CD)))
mbat_MTO = R/(eta_EM*eta_p*(e_f/g)*(CL/CD))



#fuelfractioncruise_propeller = eta_p/(g*c_p)
#fuelfractioncruise_jet = V_cruise/(g*c_j)

#W4_5 = np.exp(R*(CD/CL)*(1/fraction_propeller))     # W4/W5 for propeller
#W4_5 = np.exp(R*(CD/CL)*(1/fraction_jet))           # W4/W5 for jet

#W8/W9
#fuelfractionloiter_propeller = eta_p/(V*g*c_p)
#fuelfractionloiter_jet = 1/(g*c_j)

#W8_9 = np.exp(E*(1/fuelfractionloiter_propeller)*(CD/CL))   # W8/W9 for propeller
#W8_9 = np.exp(E*(1/fuelfractionloiter_jet)*(CD/CL))   # W8/W9 for jet


W1_TO =                      #statistics dependent on aircraft
W2_1 =                       #statistics dependent on aircraft
W3_2 =                       #statistics dependent on aircraft
W4_3 =                       #statistics dependent on aircraft
W6_5 =                       #statistics dependent on aircraft
W7_6 =                       #statistics dependent on aircraft
W8_7 =                       #statistics dependent on aircraft
W10_9 =                      #statistics dependent on aircraft
WF_10 =                      #statistics dependent on aircraft

Mff = W1_TO*W2_1*W3_2*W4_3*(1/W4_5)*W6_5*W7_6*W8_7*(1/W8_9)*W10_9*WF_10                    

Mres = 0.25
Mused = (1-Mff)
WPAX = 200*0.453592*PAX
WPAXBAGGAGE = 40*0.453592*PAX
WCargo = #Cargo Weight
WPLtot = WPAX + WPAXBAGGAGE + WCargo
MTOW = (b + WPLtot)/(Mff-a-Mres*(1-Mff))
##REWRITTEN main FORMULA & add output formulas
