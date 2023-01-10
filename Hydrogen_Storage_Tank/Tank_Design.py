import numpy as np
from Class_II_Weight_Estimation.Class_II import p_cruise
from Initial_Aircraft_Sizing.Fuselage import D_outer
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_f
from Constants import *

#Hydrogen Tank Design
#Required Volume:
LH2_vol = (m_f/71)*1.03
print('LH2 Volume', LH2_vol)
#Insulation
t_insulation = 0.14 #[m]

#Geometry
Outer_dia = D_outer
Inner_dia = Outer_dia-t_insulation
print('Inner Diameter', Inner_dia)

#part I half Sphere
v_partI = (4/3*np.pi*(Inner_dia/2)**(3))/2
s_partI = (4*np.pi*(Inner_dia/2)**(2))/2

#part III eliptical
r_flat_sphere = Inner_dia/4
ee = np.sqrt(1 - (((r_flat_sphere)**2)/((Inner_dia/2)**2)))
v_partIII = (4/3*np.pi*(Inner_dia/2)**(2)*r_flat_sphere)/2
s_partIII = (2*np.pi*(Inner_dia/2)**(2)*(1+(1-(ee)**2)/(ee*np.tanh(ee))))/2

#part II Cilinder
v_partII = LH2_vol-v_partI-v_partIII
lenght_cyl = v_partII/(np.pi*(Inner_dia/2)**2)
s_partII = 2*np.pi*(Inner_dia/2)*lenght_cyl

#Pressure
P_vent = 300000 #Pa
P_amb_fl = p_cruise
r_shell = Inner_dia/2
sigma = 172000000 #stress, Pa, r 2019-T851 aluminium alloy (ρal = 2840kg /m3) at -252◦, for 40000 cycles and with a fatigue quality index of 5
e_w = 0.8 #safetyfactor

#Accelorations
# These are Hydrostatic accelerations
rho_mean = 71 #kg/m3 #mean fuel density
K = 9 #Acceleration from loading diagram        Should probably check this value
L_tank = lenght_cyl+Inner_dia/2+r_flat_sphere #L tank in acceleration direction
Pressure_diff_acc = rho_mean*K*g*L_tank
print('pressure diff', Pressure_diff_acc)

#Structual Loads
t_shell_strucload = 0

#Tank Divider


#Results:
#Pressure
t_shell_pres = ((P_vent-P_amb_fl+Pressure_diff_acc)*r_shell)/(sigma*e_w)
print('t_shell', t_shell_pres, '[m]')

#Weight
t_tank_total = t_shell_pres + t_shell_strucload
rho_Al = 2840 #kg/m3   #2019-T851 aluminium alloy (ρal = 2840kg /m3) at -252◦
tank_surface_total = s_partI+s_partII+s_partIII
Tank_weight = tank_surface_total*rho_Al*t_tank_total
print('Tank Weight', Tank_weight)

#Geometry
print('lenght cylinder', lenght_cyl)
print('tank length', L_tank)
print('spartI', s_partI)
print('spartII', s_partII)
print('spartIII', s_partIII)
print('vpartI', v_partI)
print('vpartII', v_partII)
print('vpartIII', v_partIII)

#C.G. location starting from the forward part of the tank

