#100% fuel calculations
import numpy as np
from matplotlib import pyplot as plt

g = 9.80665

# Propulsion charecteristics
n_pk =  0.85
n_ph = 0.85
#PSFC = 0.48*(0.45/(745*3600))
e_fk = 42.9 * 10**6
e_fh = 142 * 10**6
n_engk = 0.45
n_engh = 0.3
# n_engk = (1/e_fk)*(1/PSFC)
# n_engh = (1/e_fh)*(1/PSFC)

#Aerodynamic Charecteristics
LD_crs = 16.2
L_D = 16          #Lift over drag ---> CHECK!!!

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time

# Masses Aircraft
m_mto = 19314   # in kg --> Max take off
m_oe =  12486   # in kg --> Operating empty
m_pldes = 5730  # in kg --> Design payload
m_plmax = m_pldes*1.1   # in kg --> Max payload
#m_f = 1098 # in kg --> Fuel mass
# #Point A
ranges1 = [0]
plmasses1 = [m_plmax]

#Point B (Range at Max Payload)

#m_fB = m_mto * (1 - np.exp(-R / (n_engh * n_ph * (e_fh /g) * (L_D))))
m_fB = m_mto - m_oe - m_plmax

R_b1 = n_engh * n_ph * (L_D) * (e_fh/g) * np.log((m_oe + m_plmax + m_fB)/(m_oe + m_plmax))
R_lostB = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eqB = ((R_b1 + R_lostB)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
R_auxB = R_eqB - R_b1
R_b = R_b1 - R_auxB

ranges1 = np.append(ranges1, [R_b])
plmasses1 = np.append(plmasses1, [m_plmax])

#print(ranges1)
#print(R_auxB)
#print(R_eqB)
#print(R_lostB)

#Point C (Design Range)
R_lost1 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eq1 = ((R_nom + R_lost1)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
R_aux1 = R_eq1 - R_nom
#m_fC = m_mto * (1 - np.exp(-R_nom / (n_engh * n_ph * (e_fh /g) * (L_D))))
m_fC = m_mto - m_pldes - m_oe
R_c1 = (n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_fC)/(m_oe + m_pldes))) - R_aux1

ranges1 = np.append(ranges1, [R_c1])
plmasses1 = np.append(plmasses1, [m_pldes])

#print("Rc", R_c1)
print("mf" ,m_fC)
#print(R_eq1)
#print(R_aux1)
#print(R_lost1)

#Point D (Ferry Range)
R_d1 = (n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe  + m_fC)/(m_oe ))) - R_aux1

ranges1 = np.append(ranges1, [R_d1])
plmasses1 = np.append(plmasses1, [0])

#print(R_d1)

#Constructing the actual plot

# plotting the points
plt.plot(ranges1, plmasses1, color='green', linewidth=2,
         marker='o', markerfacecolor='green', markersize=5)


# print(ranges1)
# print(plmasses1)
# print(points)
# naming the x axis
plt.xlabel('Range (m)')
# naming the y axis
plt.ylabel('Payload (Kg)')

# giving a title to my graph
plt.title('Payload-Range diagram for LH2 combustion')

# function to show the plot
plt.show()

#------------

#Hydrogen- fuel case


# # #Point A
# ranges2 = [0]
# plmasses2 = [m_pldes]
# #Point B (Range at Max Payload)
# # R_b = n_eng * n_p * (L/D) * (e_f /g) * np.log((m_oe + m_plmax + m_f)/(m_oe + m_pl))
# #
# # m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))
# #
#
# #m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))
#
# ranges = np.append(ranges, [R_b])
# plmasses = np.append(plmasses, [m_pl])
#
# print(ranges)
#
# #Point C (Design Range)
# ratio_h = np.linspace(0.0, 1.0, num=21)
# ratio_k = 1 - ratio_h
# f_con = 5/100
# R_div = 200000
# t_E = 45 * 60
# R_lost2 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
# print(R_lost2)
# R_nom2 = 1852000
# R_eq2 = ((R_nom2 + R_lost2)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
# print(R_eq2)
# R_aux2 = R_eq2 - R_nom2
# print(R_aux2)
# m_mto = 21500
# m_fh = ((m_mto * (1 - np.exp(-R_nom2 / (n_engh * n_ph * (e_fh /g) * (L_D)))))*0.3)
# m_fk =  ((m_mto * (1 - np.exp(-R_nom2 / (n_engk * n_pk * (e_fk /g) * (L_D)))))*0.7)
# R_c2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_fh)/(m_oe + m_pldes))) +  (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe + m_pldes + m_fk)/(m_oe + m_pldes))))   - R_aux1
# print("Rc", R_c2)
# ranges2 = np.append(ranges2, [R_c2])
# plmasses2 = np.append(plmasses2, [m_pldes])
#
# #Point D (Ferry Range)
#
# R_d2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe  + m_fh)/(m_oe )))+ (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe  + m_fk)/(m_oe )))) - R_aux1
# print(R_d2)
#
# ranges2 = np.append(ranges2, [R_d2])
# plmasses2 = np.append(plmasses2, [0])
#
# #Constructing the actual plot
#
# # plotting the points
# plt.plot(ranges2, plmasses2, color='red', linewidth=2,
#          marker='o', markerfacecolor='red', markersize=5)
#
#
# # naming the x axis
# plt.xlabel('Range (m)')
# # naming the y axis
# plt.ylabel('Payload (Kg)')
#
# # giving a title to my graph
# plt.title('Payload-Range diagram for LH-Kerosene hybrid')
#
# # function to show the plot
# plt.show()
#
# a = np.linspace(0.0, 1.0, num=21)
# print(a)
# b = 1-a
#
# print(b)
