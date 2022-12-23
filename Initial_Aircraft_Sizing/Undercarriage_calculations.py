import numpy as np
from matplotlib import pyplot as plt
def myround(x, base=5):
    return base * round(x/base)

#Variables
#Wheel tyres
N_nw = 2
W_TO = 23000*9.81
LCN = 19                            #value used for Fokker 27

#Clearance Angles
overturn = 55* np.pi *(1/180)


span = 27.06
y_e = 4.05     #distance to outmost engine, (span of wing for LHFC, EHSP, EHS and 4.05m for H2C
z_t = 4.0108     #wingtip height  (4.0108m if HW, 1.0108 if LW)
z_n = 1.015    #outmost engine clearance
phi = 5* np.pi *(1/180)

#Disposition Calculation Variables
l_fus = 26.876
x_cg = 12.04
l_tail = 9.03
X_tailtocg =  l_fus - x_cg - l_tail
D_finner = 3.01
scrape = 18* np.pi *(1/180)
tipback = 18* np.pi *(1/180)

#Calculation of Nr of Landing gears needed
N_mw = W_TO / 60000
N_mw = myround(N_mw, base=4)

print("Number of main wheels", N_mw)


#Calculation of tyre pressure
p = 430*np.log(LCN) - 680
print("Infl pressure", p, (p*0.0101972))

#Wheel sizing: Dimensions

P_mw = (0.92*W_TO) / N_mw       #Static load per main gear wheel
P_nw = (0.08*W_TO) / N_nw       #Static load per nose gear wheel

print("P_mw =", P_mw, P_mw/9.81)
print("P_nw =", P_nw, P_nw/9.81)
# #Now tire is selected from Torenbeek's plots
#
#
#
# #Undercarriage disposition PAPER DRAWING AID RECOMMENDED
x_aftcg = x_cg
a = X_tailtocg * np.sin(18*np.pi /180)
b = D_finner / 2
l_mlg = (a+b)* np.cos(18*np.pi /180)* np.sin(18*np.pi /180)
z_mlg = (a+b)*np.cos(18*np.pi /180)*np.cos(18*np.pi /180)
x_mlg = l_mlg + x_aftcg
h_mlg = z_mlg - b
print(X_tailtocg)
print(a)
print("lmlg", l_mlg)
print("xmlg", x_mlg)
print("zmlg", z_mlg)
print("hmlg", h_mlg)

#Nosewheel placement

l_nlg = (0.92/0.08)*l_mlg
print("lnlg", l_nlg)
x_nlg = x_aftcg - l_nlg

print("xnlg", x_nlg)


#Dummy values
#Clearances
#Lateral tip-over clearance


#Lateral tip-over clearance y_MLG > y_MLG_tov
y_MLG_tov = (l_nlg + l_mlg)/ np.sqrt(((l_nlg **2 * (np.tan(overturn))**2)/(z_mlg)**2)-1)
print("yMLGtov", y_MLG_tov)


#Tip Clearance  y_MLG > y_MLG_tip

y_MLG_tip = span/2 - z_t / np.tan(phi)
print("yMLGtip", y_MLG_tip)

#Engine Clearance  y_MLG > y_MLG_eng

y_MLG_eng = y_e - z_n / np.tan(phi)

print("yMLGeng", y_MLG_eng)
