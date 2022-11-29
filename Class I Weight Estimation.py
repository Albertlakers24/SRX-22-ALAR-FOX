import numpy as np
# Outputs: MTOW OEW and Wf
# Input: Range, Cruise velocity,
# Relationships: MTOW = WOE + WF + WPL
# Relationship: WOE = WE + Wtfo + Wcrew
# Relationship: Wtfo = Mtfo * Wto
# Relationship: WOE = a * MTOW + b
# Cd = 2*Cd0
# Cl = np.sqrt(np.pi()*A*Cd0 * A)
#
a = input("What is the a value of the linear regression?")
b = input("What is the b value of the linear regression?")
g = 9.80665
A = 10#Aspect ratio (CHANGE)
Cd0 = 2 #zero lift drag (CHANGE)
R = 1000 * 1852 # range in meters
E = 35 * 60 # loiter endurance in seconds
