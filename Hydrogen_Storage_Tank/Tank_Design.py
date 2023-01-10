import numpy as np
from Class_II_Weight_Estimation.Class_II import p_cruise, W_ver_tail, W_hor_tail, W_fus
from Initial_Aircraft_Sizing.Fuselage import D_outer, l_t, l_f
from Initial_Aircraft_Sizing.Empennage_Design import x_h, x_v
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_f
from Constants import *
rho_Al = 2840 #kg/m3   #2019-T851 aluminium alloy (ρal = 2840kg /m3) at -252◦

#Hydrogen Tank Design
#Required Volume:
LH2_vol = (m_f/71)*1.03
print('LH2 Volume', LH2_vol)

#Insulation
t_insulation = 0.1527 #[m] #Should check this from the Boil-off and thus mass flow

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
print('pressure diff by acc', Pressure_diff_acc)

#Structual Loads
#print('lf', l_f)
dis_tank_ver =l_t-r_flat_sphere-(l_f-x_v)
dis_tank_hor =l_t-r_flat_sphere-(l_f-x_h)
load_bending = g*(W_hor_tail*dis_tank_hor+W_ver_tail*dis_tank_ver+W_fus/l_f*l_t*l_t/2)
t_bending = (load_bending/sigma)/(np.pi*Inner_dia)
t_shell_strucload = t_bending*4 #Loading diagram max acc
#print('t_bend',t_shell_strucload)

#Tank Divider
# Assumed to be the same as the forward cap #TU delft report
# Added in weight part

#Results:
'''
Things to check: Insulation thickness
Accelorations
Structual Loads
'''
#Pressure
t_shell_pres = ((P_vent-P_amb_fl+Pressure_diff_acc)*r_shell)/(sigma*e_w)
#print('t_shell', t_shell_pres, '[m]')

#Weight
t_tank_total_SF = (t_shell_pres + t_shell_strucload)*1.5 #With Safety Factor
print('t_total', t_tank_total_SF)
tank_surface_total = s_partI+s_partII+s_partIII
tank_divider_weight = s_partIII*rho_Al*t_tank_total_SF #Divider
Tank_weight = tank_surface_total*rho_Al*t_tank_total_SF+tank_divider_weight
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
cg_tank = ((s_partIII*r_flat_sphere/2)+(s_partII*(r_flat_sphere+lenght_cyl/2))+(s_partI*(r_flat_sphere+lenght_cyl+Inner_dia/4))+(s_partIII*1.5*r_flat_sphere))/(tank_surface_total+s_partIII)
print('C.G. Tank', cg_tank)

#Boil-Off
k_mli = 1.43 #W/m2   #MLI Heat transfer rate                #0.135 * 0.001 #W / m * K
delta_h = 461000 #J/kg
mass_flow = 0.12 #kg/s
Heat_flow_allowed = mass_flow * delta_h  #W #Input
T_ex_max = 323 #K
T_in_min = 20 #K
Insulation_resistance = (T_ex_max-T_in_min)/Heat_flow_allowed
Insulation_thickness = Insulation_resistance*k_mli*(s_partIII+s_partII+s_partI)
print('Tank Insulation Thickness', Insulation_thickness)

