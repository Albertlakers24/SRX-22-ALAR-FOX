import numpy as np

# Outputs: MTOW OEW and WF
# Relationships: MTOW = WOE + WF + WPL
# Relationships: MTOW = WE + WF + WPLtot + Wtfo
# LOOK INTO MF/MTO FORMULA FROM ROLOEF READER
# LIQUID FUELS WILL HAVE THE SAME FORMULA AS JET-A FUEL
# BATTERIES WILL HAVE FORMULA WITHOUT ln()
#

#Constants
g = 9.80665
R = 1000 * 1852                 #Range in meters
E = 35 * 60                     #Loiter endurance in seconds
Mtfo = #0.001 - 0.005           #Trapped fuel oil in fraction
Wcrew= 3*190*0.45359237*g       #in N
V_cruise =                      #m/s


#Constant on estimation
A =                             #Aspect ratio (CHANGE) - high A for slender wing
e =                             #Oswald factor (CHANGE)
WPL = *g                        #From the guidelines
eta_p =                         #propeller efficiency   -> maximize
c_p =                           #propeller              -> minimize
c_j =                           #jet                    -> minimize

WPLtot = WPL + Wcrew

#linear regression relating to MTOW depending on the aircraft   WE = a*MTOW+b
a = 0.5422 or 0.4985            #Linear regression from MTOW/ OEW (turboprop or turbojet)
b = 1455.2 or 1782.3            #Linear regression from MTOW/ OEW (turboprop or turbojet)

#Wtfo = Mtfo * MTOW    (trapped fuel and oil)


#WF use the fuel fraction method
#WF = Mused*MTOW*(Mres + 1)

#Cdo calculations
Cfe =                           #Equivalent skin friction coefficient - depending on aircraft from empirical estimation
Swet_S =                        #Wetted area ratios - depending on airframe structure
Cd0 = Cfe*Swet_S

#W4/W5
CL = np.sqrt(np.pi()*Cd0*A*e)
CD = 2 * Cd0

fuelfractioncruise_propeller = eta_p/(g*c_p)
fuelfractioncruise_jet = V_cruise/(g*c_j)

W4_5 = np.exp(R*(CD/CL)*(1/fraction_propeller))     # W4/W5 for propeller
W4_5 = np.exp(R*(CD/CL)*(1/fraction_jet))           # W4/W5 for jet

#W8/W9
fuelfractionloiter_propeller = eta_p/(V*g*c_p)      
fuelfractionloiter_jet = 1/(g*c_j)

W8_9 = np.exp(E*(1/fuelfractionloiter_propeller)*(CD/CL))   # W8/W9 for propeller
W8_9 = np.exp(E*(1/fuelfractionloiter_jet)*(CD/CL))   # W8/W9 for jet

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


##REWRITTEN main FORMULA & add output formulas
