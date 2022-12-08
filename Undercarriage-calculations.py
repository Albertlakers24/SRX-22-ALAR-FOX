import numpy as np
from matplotlib import pyplot as plt
def myround(x, base=5):
    return base * round(x/base)

#Variables
#Wheel tyres
N_nw = 2
W_TO = 180000
LCN = 19    #value used for Fokker 27

#Clearance Angles
overturn = 33
l_n = 2
l_m = 2
z = 1
b = 5
y_e = 3
z_t = 2
z_n = 1
phi = 11


#Calculation of Nr of Landing gears needed
N_mw = W_TO / 60000
N_mw = myround(N_mw, base=4)

print(N_mw)


#Calculation of tyre pressure
p = 430*np.log(LCN) - 680
print(p)

#Wheel sizing: Dimensions

P_mw = (0.92*W_TO) / N_mw       #Static load per main gear wheel
P_nw = (0.08*W_TO) / N_nw       #Static load per nose gear wheel

print("P_mw =", P_mw)
print("P_nw =", P_nw)
#Now tire is selected from Torenbeek's plots



#Undercarriage disposition should be done on paper.]






#Dummy values
#Clearances
#Lateral tip-over clearance


#Lateral tip-over clearance y_MLG > y_MLG_tov
y_MLG_tov = (l_n + l_m)/ np.sqrt(((l_n **2 * (np.tan(overturn))**2)/z**2)-1)



#Tip Clearance  y_MLG > y_MLG_tip

y_MLG_tip = b/2 - z_t / np.tan(phi)


#Engine Clearance  y_MLG > y_MLG_eng

y_MLG_eng = y_e - z_n / np.tan(phi)