from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *

#all the equations from Raymer

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
tc = 0.12

#Ulitimate load = 1.5*limit load factor (assumption for now)
Lim_load = 4.
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
W_hor_tail = W_hor_tail_lbs *lbs_kg

#Vertical Tail
#Ht/Hv
H_t_H_v = 1.0
#S_vt
S_vt = Sv/(ft_m)**2
#taper_vt
taper_vt = taperv

W_ver_tail_lbs = 0.073*(1+0.2*(H_t_H_v))*(N_z*W_dg)**(0.376)*q_lbft2**(0.122)*S_vt**(0.873)*((100*tc)/np.cos(Lambda))**(-0.49)*(A/(np.cos(Lambda))**(2))**(0.357)*taper_vt**(0.039)
W_ver_tail = W_ver_tail_lbs *lbs_kg


#Fuselage
from Initial_Aircraft_Sizing.Fuselage import S_f_wet
from Class_I_Weight_Estimation.Drag_Polar import max_CL_CD

#S_f
S_f = S_f_wet/(ft_m)**2

#L_t
L_t = x_h/ft_m

#L_D
L_D = max_CL_CD

#W_press
V_pr_m3 = 100 #not right!!!!
V_pr_ft3 = V_pr_m3/ft_m**(3)
P_delta = 8 #psi
W_press = 11.9+(V_pr_ft3*P_delta)**(0.271)

W_fus_lbs = 0.052*S_f**(1.086)*(N_z*W_dg)**(0.177)*L_t**(-0.051)*(L_D)**(-0.072)*q_lbft2**(0.241)+W_press
W_fus = W_fus_lbs * lbs_kg

#Main Landing Gear
#N_l load factor
N_main_load = 4 #not right!!
N_l = 1.5* N_main_load

#W_l Max land weight,lbs
W_l =  beta_s_land_fc*m_mto/lbs_kg
#L_m lenght of main landing gear
L_m = 1.2/ft_m #not right!!


W_land_main_lbs = 0.095*(N_l*W_l)**0.768*(L_m/12)**(0.409)
W_land_main = W_land_main_lbs * lbs_kg

#Nose Landing Gear

#Installed Engine
#(From our own estimates and NOT Raymer)

#Fuel System
#(From our own estimates and NOT Raymer)

#Flight Controls

#Hydrolics

#Electrical

#Avionics

#Airconditioning and Anti-Icing

#Furnishing

#Results
print('Weight Wing kg =', W_wing)
print('Weight Horizontal Tail = ', W_hor_tail)
print('Weight Vertical Tail = ', W_ver_tail)
print('Weight Fuselage =', W_fus)
print('Weight Main Landing Gear', W_land_main)


