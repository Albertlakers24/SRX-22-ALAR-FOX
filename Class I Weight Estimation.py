# Output is MTOW and OEW
#Outputs: WTO, WOE, WF, WFused

# Input is
# Relationships: MTOW = WOE + WF + WPL
# Relationship: WOE = WE + Wtfo + Wcrew

##FUEL Calculations

# Relationship:
# Wtfo = Mtfo * MTOW    (trapped fuel and oil) -> Mtfo 0.001 - 0.005
# WPLtot = WPL + Wcrew
# MTOW = WE + WF + WPLtot + Mtfo*WTO
# WF = WFres + WFused
# WFres =
# WFused = Mused + MTOW
# Mused = 1-Mff
# Mff is the fraction method
# R = (/)*(L/D)*ln(W4/W5) dependent on fueltype/propulsion system
# E = (/)*(L/D)*ln(W8/W9)

# FUEL DEPENDENT:
# * Mtfo
# * Mff


#Inputs:    Dependent on fuel-type
#           statics on Mff to be found on the slides (except cruise and loiter):
#           for W4/WF5 (cruise) -> R      --to be calculated
#           for W8/W9 (loiter) -> E
#           Payload from statistics
#           MTOW = a*WF +b ---given by Gabriel
#           Mtfo = ?
#           Wcrew = ?

#WOE
# WOE = WE + Wtfo + Wcrew
#
# Wtfo = Mtfo * MTOW    (trapped fuel and oil) -> Mtfo 0.001 - 0.005
g = 9.80665
Mtfo = #0.001 - 0.005
Wcrew= 3*190*0.45359237*g     #in N
WE = ##linear regression relating to MTOW depending on the aircraft


##WF
#WF = WFres + WFused
#WFres = Mres *WFused
#WFused = Mused*MTOW


#WF = Mused*MTOW*(Mres + 1)

Mres = 0.25
Mff = ##from the fuel fraction method
Mused = (1-Mff)








