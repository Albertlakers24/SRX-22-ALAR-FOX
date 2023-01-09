

from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto, m_tanks, oem
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *
from Constants import *

T_cruise, p_cruise, rho_cruise, a_cruise = ISA_calculator(h_cruise, 0)

#all the equations from Raymer
#Calculated seperaterly:

#LH2 Storage
#Estimated at 1.4 grav index, add more detail later I would say
LH2_system_tank = m_tanks             #Must be done, thus change

#Engine
power_req = m_mto * g / W_P_design
Eng_W_kg = 6000 #Check Later
Engine_weight = power_req / Eng_W_kg

#Fuel Cell
FC_W_kg = 3000
Fuel_Cell_Weight = power_req / FC_W_kg

#Wing
#Surface, ft^2
S = m_mto * g / W_S_design
S_w = S/((ft_m)**2)

#weight of fuel in wing, lb = 0
W_fw = 1.
#aspect ratio
A = 12.
#Sweep = 0
Lambda = 0
#dynamic pressure at cruise lb/ft^2
rho_cruise = ISA_calculator(280*FL_ft*ft_m,dt_cruise)[2]
q_kgm2 = 0.5*rho_cruise*V_cruise**(2)
q_lbft2 = q_kgm2*23.73 #kg/m2 to lbs/ft2

#taper ratio
taper = 0.45

#t/c (assumption for now)
tc = 0.12           #Change maybe??

#Ulitimate load = 1.5*limit load factor (assumption for now)
Lim_load = 4.       #Will change, add link to maneuver load diagram
N_z = 1.5*Lim_load

#Design gross weight, lb (MTOW)
W_dg = m_mto/lbs_kg


W_wing_lbs = 0.036*S_w**(0.758)*W_fw**(0.0035)*(A/(np.cos(Lambda))**(2))**(0.6)*q_lbft2**(0.006)*taper**(0.04)*((100*tc)/np.cos(Lambda))**(-0.3)*(N_z*W_dg)**(0.49)
W_wing = W_wing_lbs * lbs_kg


#Horizontal Tail
from Initial_Aircraft_Sizing.Empennage_Design import *
#S_ht
S_ht = Sh/(ft_m)**2
#taper_ht
taper_ht = taperh

W_hor_tail_lbs = 0.016*(N_z*W_dg)**(0.414)*q_lbft2**(0.168)*S_ht*(0.896)*((100*tc)/np.cos(Lambda))**(-0.12)*(A/(np.cos(Lambda))**(2))**(0.043)*taper_ht**(-0.02)
W_hor_tail = W_hor_tail_lbs * lbs_kg

#Vertical Tail
#Ht/Hv
H_t_H_v = 1.
#S_vt
S_vt = Sv/(ft_m)**2
#taper_vt
taper_vt = taperv

W_ver_tail_lbs = 0.073*(1+0.2*(H_t_H_v))*(N_z*W_dg)**(0.376)*q_lbft2**(0.122)*S_vt**(0.873)*((100*tc)/np.cos(Lambda))**(-0.49)*(A/(np.cos(Lambda))**(2))**(0.357)*taper_vt**(0.039)
W_ver_tail = W_ver_tail_lbs *lbs_kg


#Fuselage
from Initial_Aircraft_Sizing.Fuselage import S_f_wet, l_f
from Class_I_Weight_Estimation.Drag_Polar import max_CL_CD

#S_f
S_f = S_f_wet/(ft_m)**2

#L_t
L_t = x_h/ft_m

#L_D
L_D = max_CL_CD

#W_press
V_pr_m3 = 86.78 #Could change
V_pr_ft3 = V_pr_m3/(ft_m**(3))
P_delta = 8 #psi
W_press = 11.9+(V_pr_ft3*P_delta)**(0.271)

W_fus_lbs = 0.052*S_f**(1.086)*(N_z*W_dg)**(0.177)*L_t**(-0.051)*(L_D)**(-0.072)*q_lbft2**(0.241)+W_press
W_fus = W_fus_lbs * lbs_kg

#Main Landing Gear

#W_l Max land weight,lbs
W_l =  beta_s_land_fc*m_mto/lbs_kg
#L_m lenght of main landing gear, in
L_m = 0.853*12/ft_m
#N_l load factor
N_main_load = 240000/(g*beta_s_land_fc*m_mto) #240000 is from 4 * 60000
N_l = 1.5* N_main_load

W_land_main_lbs = 0.095*(N_l*W_l)**(0.768)*(L_m/12)**(0.409)
W_land_main = W_land_main_lbs * lbs_kg

#Nose Landing Gear
#L_n lenght of nose landing gear, in
L_n = 0.853*12/ft_m

W_land_nose_lbs = 0.125*(N_l*W_l)**(0.566)*(L_n/12)**(0.845)
W_land_nose = W_land_nose_lbs * lbs_kg

#Installed Engine
N_engines = 1           #Could change
W_installed_eng = 2.575*Engine_weight**(0.922)*N_engines

#Fuel System
#this is done with the gravimetric index etc and hydrogen tank, Should be a better approach

#Flight Controls
#L_fus
L_fus = l_f/ft_m
#B_w wing span
B_w = b/ft_m

W_flight_controls_lbs = 0.053*L_fus**(1.536)*B_w**(0.371)*(N_z*W_dg*0.0001)**(0.80)
W_flight_controls = W_flight_controls_lbs * lbs_kg


#Hydrolics
W_hydraulics_lbs = 0.001*W_dg
W_hydraulics = W_hydraulics_lbs * lbs_kg

#Avionics
#W_auv uninstalled avionics weight
W_uav = 800

W_av_lbs = 2.117*W_uav**(0.933)
W_av = W_av_lbs * lbs_kg

#Electrical
W_elec_lbs = 12.57*(W_av_lbs)**(0.51)
W_elec = W_elec_lbs * lbs_kg

#Airconditioning and Anti-Icing
#N_p personnel
N_p = 3
#Mach Number
Mach = V_cruise/a_cruise

W_aircon_ice_lbs = 0.265*W_dg**(0.52)*N_p**(0.68)*W_av_lbs**(0.17)*Mach**(0.08)
W_aircon_ice = W_aircon_ice_lbs * lbs_kg

#Furnishing
W_furnish_lbs = 0.0582*W_dg-65
W_furnish = W_furnish_lbs * lbs_kg



#Results:

#Things to be reconsidered later before iteration etc:
#LH2 storgare mass
#Engine mass
#t/c ratio (both wing and empennage)
#Lim load
#Fuselage pressurized volume
#Landing gear load
#Both landing gear length

print('Results Start Here:')
print('Weight Storage System LH2', LH2_system_tank)
print('Weight Engines =', Engine_weight)
print('Weight Fuel Cell', Fuel_Cell_Weight)
print('Weight Wing =', W_wing)
print('Weight Horizontal Tail = ', W_hor_tail)
print('Weight Vertical Tail = ', W_ver_tail)
print('Weight Fuselage =', W_fus)
print('Weight Main Landing Gear =', W_land_main)
print('Weight Nose Landing Gear =', W_land_nose)
print('Weight Flight Controls =', W_flight_controls)
print('Weight Hydraulics =', W_hydraulics)
print('Weight Avionics =', W_av)
print('Weight Electronics =', W_elec)
print('Weight Aircon and De-Icing', W_aircon_ice)
print('Weight Furnishing =', W_furnish)

print('Class II Weight Estimation =', LH2_system_tank+Engine_weight+Fuel_Cell_Weight+W_wing+W_hor_tail+W_ver_tail+W_fus+W_land_main+W_land_nose+W_flight_controls+W_hydraulics+W_av+W_elec+W_aircon_ice+W_furnish)
print('Class I OEM =', oem)



