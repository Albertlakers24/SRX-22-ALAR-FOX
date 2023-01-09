import numpy as np
from Class_II_Weight_Estimation.Class_II import p_cruise
from Initial_Aircraft_Sizing.Fuselage import D_outer
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_f

#Hydrogen Tank Design
#Required Volume:
LH2_vol = (m_f/71)*1.03

#Insulation
t_insulation = 0.14

#Geometry
Outer_dia = D_outer
Inner_dia = Outer_dia-t_insulation
print(Inner_dia)
#part I half Sphere
v_partI = (4/3*np.pi*(Inner_dia/2)**(3))/2
s_partI = (4*np.pi*(Inner_dia/2)**(2))/2

#part II Cilinder

#part III eliptical
r_flat_sphere = Inner_dia/4
ee = np.sqrt(1 - (((r_flat_sphere)**2)/((Inner_dia/2)**2)))
v_partIII = (4/3*np.pi*(Inner_dia/2)**(2)*r_flat_sphere)/2
s_partIII = 2*np.pi*(Inner_dia/2)**(2)*(1+(1-(ee)**2))/ee*np.tanh(ee)**(-1)/2

print('spartI', s_partI)
print('spartIII', s_partIII)
print('vpartI', v_partI)
print('vpartIII', v_partIII)
#Pressure
P_vent = 300000 #Pa
P_amb_fl = p_cruise
r_shell = Inner_dia/2
sigma = 172000000 #stress, Pa, r 2019-T851 aluminium alloy (ρal = 2840kg /m3) at -252◦, for 40000 cycles and with a fatigue quality index of 5
e_w = 0.8 #safetyfactor
t_shell = ((P_vent-P_amb_fl)*r_shell)/(sigma*e_w)
print('t_shell', t_shell, '[m]')

#Accelorations

#Structual Loads

#Tank Divider


#Results:
#Weight

#Geometry
