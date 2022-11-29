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
a = 0.5422 or 0.4985 #Linear regression from MTOW/ OEW (turboprop or turbrojet)
b = 1455.2 or 1782.3 #Linear regression from MTOW/ OEW (turboprop or turbrojet)
g = 9.80665
A = #Aspect ratio (CHANGE)
Cd0 = #zero lift drag (CHANGE)
R = 1000 * 1852 # range in meters
E = 35 * 60 # loiter endurance in seconds
