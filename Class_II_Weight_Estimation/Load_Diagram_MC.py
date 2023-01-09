import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto, oem, m_zf, m_f
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *
from Initial_Aircraft_Sizing.Wing_planform import b, c_mac

# Density
rho = 0 # just to define rho, this NEVER changes
density = 0 # 0 @ sea-level, 1 @ cruise, else @ loiter
if density == 0:
    rho = rho_0 # density at sea level
elif density == 1:
    rho = rho_cruise
else:
    rho = 0.663838

#MANEUVER DIAGRAM DESIGN
#Max lift coefficient
if density == 0:
    Clmax = CL_max_landing
elif density == 1:
    Clmax = CL_max_cruise
else:
    Clmax = CL_max_loiter # add loiter Clmax to Constants file

#Load factor
m_mto_lbs = m_mto / lbs_kg #conversion between kg and lbs
n = 2.1 + (24000 / (m_mto_lbs + 10000))

if n < 2.5:
    n_max = 2.5
elif n > 3.8:
    n_max = 3.8
else:
    n_max = n

n_min = -1

#Cruise speed EAS
V_C = V_cruise
V_C = V_C * np.sqrt(rho/rho_0)

#Stall speed EAS
V_S = np.sqrt((2 * 1 * W_S_design) / (rho * Clmax)) * np.sqrt(rho / rho_0)

#Dive speed EAS
if density == 0:
    h = 0
    dt = dt_takeoff
elif density == 1:
    h = h_cruise
    dt = dt_cruise
else:
    h = h_loiter
    dt = dt_loiter

a = ISA_calculator(h,dt)[3]
M_C = V_C / a
M_D = M_C + 0.05

V_D1 = V_C/0.8
V_D2 = M_D * a
V_D = min(V_D1, V_D2)

#---- CONSTRUCTION OF MANEUVER DIAGRAM ----
#From 0 to V_A
V_A = []
n_A_list = []
V_A_fin = V_S * np.sqrt(n_max)
for i in np.arange(0, V_A_fin, 0.1):
    V_A.append(i)
    n = 0.5 * i**2 * rho_0 * Clmax / W_S_design
    n_A_list.append(n)

#From 0 to V_H
V_H = []
n_H_list = []
for j in np.arange(0, V_S, 0.1):
    V_H.append(j)
    n = -0.5 * j**2 * rho_0 * Clmax / W_S_design
    n_H_list.append(n)
'''
#For flaps
n_max_flaps = 2
n_flaps_list = []
V_flaps = []
V_max_flaps = np.sqrt((2 * n_max_flaps * W_S_design) / (rho_0 * CL_max_landing)) #double check if constant rho_0 or not
V_intersect = np.sqrt((2 * n_max_flaps * W_S_design) / (rho_0 * CL_max_cruise)) #double check if constant rho_0 or not
for i in np.arange(0, V_max_flaps, 0.1):
    V_flaps.append(i)
    n = (1/2 * CL_max_landing * rho_0) * i**2 / W_S_design #double check if constant rho_0 or not
    n_flaps_list.append(n)
'''

#Graph
point_F = [V_C, -1]
point_E = [V_D, 0]
point_A = [V_A_fin, n_max]
point_H = [V_S, -1]
x_values_FE = [point_F[0], point_E[0]]
y_values_FE = [point_F[1], point_E[1]]
point_D = [V_D, n_max]
x_values_DE = [point_D[0], point_E[0]]
y_values_DE = [point_D[1], point_E[1]]
x_values_AD = [point_A[0], point_D[0]]
y_values_AD = [point_A[1], point_D[1]]

fig, ax = plt.subplots()
plt.tick_params(labelbottom = False)
plt.plot(V_A, n_A_list, "black")
plt.plot(V_H, n_H_list, "black")
#plt.plot(V_flaps, n_flaps_list, "black")
#plt.plot([V_max_flaps, V_intersect], [2, 2], "black", marker="o")
plt.plot(x_values_AD, y_values_AD, "black", marker="o")
plt.plot(x_values_FE, y_values_FE, 'black', marker="o")
plt.plot(x_values_DE, y_values_DE, "black", marker="o")
plt.plot([point_H[0], point_F[0]], [point_H[1], point_F[1]], "black", marker="o")
plt.plot([0, V_D], [1, 1], 'black', linestyle="--")
plt.plot([V_A_fin, point_A[0]], [0, point_A[1]], 'black', linestyle="--")
plt.plot([V_S, V_S], [0, 1], "black", linestyle="--")
plt.plot([point_F[0], V_C], [point_F[1], 0], "black", linestyle="--")
plt.plot([0,0],[0,0],"black",marker="o")
plt.text(0-0.015, 0+0.25, "")
plt.text(point_A[0]-2, point_A[1]+0.2, "A")
plt.text(point_D[0]-2, point_D[1]+0.25, "D")
plt.text(point_E[0]+2, point_E[1]+0.25, "E")
plt.text(point_F[0]-1, point_F[1]-0.25, "F")
plt.text(point_H[0]-2, point_H[1]-0.25, "H")
plt.text(V_S-2,0-0.2, r'$V_{S}$')
plt.text(V_A_fin-2,0-0.2, r'$V_{A}$')
plt.text(V_C-2,0+0.1, r'$V_{C}$')
plt.text(V_D-2,0-0.25, r'$V_{D}$')
plt.xlim(left=0, right=(V_D+20))
plt.ylim(top=3,bottom=-1.5)
# plt.xlabel("V")
ax.set_xlabel("V", weight = "bold")
ax.xaxis.set_label_coords(0.96, 0.3)
plt.ylabel("n", weight = "bold")
ax.spines['bottom'].set_position(('data',0))
plt.grid()
plt.show()

#GUST DIAGRAM DESIGN
#Constants
CL_alpha = 5.03
chord = b / A
S = m_mto * g / W_S_design
W_S_oem = oem * g / S

#Gust Profile
if density == 0:
    U_ref_C = 17.07
    U_ref_D = U_ref_C / 2
elif density == 1:
    U_ref_C = 11.37
    U_ref_D = U_ref_C / 2
else:
    U_ref_C = 12.71
    U_ref_D = U_ref_C / 2

Z_mo = h_cruise
R1 = beta_s_land_fc
R2 = m_zf/m_mto
F_gz = 1 - (Z_mo / 76200)
F_gm = np.sqrt(R2 * np.tan((np.pi * R1)/4))
F_g = 0.5 * (F_gz + F_gm)
H1 = c_mac * 12.5
H2 = 107
if H1 < H2:
    H = H2
else:
    H = H1