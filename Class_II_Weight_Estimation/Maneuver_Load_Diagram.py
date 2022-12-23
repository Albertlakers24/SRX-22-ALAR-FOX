import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *
from Initial_Aircraft_Sizing.Wing_planform import b

rho = rho_0 # IF CALCULATING @ SEA-LEVEL, rho = rho_0. IF CALCULATING @ CRUISE, rho = rho_cruise

# T_cruise, p_cruise, rho_cruise, a_cruise = ISA_calculator(5000*ft_m, 0)
# V_cruise = V_cruise * np.sqrt(rho_cruise / rho_0)
constant_q = 1/2 * CL_max_cruise * rho
mass = m_mto
S = m_mto * g / W_S_design
V_A = []
n_A_list = []
n_max = max(2.1 * (24000 / (m_mto + 10000)), 2.5)
n_min = -1
V_max = np.sqrt((2 * n_max * W_S_design) / (rho_0 * CL_max_cruise)) * np.sqrt(rho / rho_0)
for i in np.arange(0, V_max, 0.1):
    V_A.append(i)
    n = constant_q * i**2 / W_S_design
    n_A_list.append(n)

V_H = []
n_H_list = []
V_min = np.sqrt((2 * 1 * W_S_design) / (rho_0 * CL_max_cruise)) * np.sqrt(rho / rho_0)
for j in np.arange(0, V_min, 0.1):
    V_H.append(j)
    n = -constant_q * j ** 2 / W_S_design
    n_H_list.append(n)

V_D = V_cruise / 0.8

n_max_flaps = 2
n_flaps_list = []
V_flaps = []
V_max_flaps = np.sqrt((2 * n_max_flaps * W_S_design) / (rho_0 * CL_max_landing)) #double check if constant rho_0 or not
V_intersect = np.sqrt((2 * n_max_flaps * W_S_design) / (rho_0 * CL_max_cruise)) #double check if constant rho_0 or not
for i in np.arange(0, V_max_flaps, 0.1):
    V_flaps.append(i)
    n = (1/2 * CL_max_landing * rho_0) * i**2 / W_S_design #double check if constant rho_0 or not
    n_flaps_list.append(n)

point_F = [V_cruise, -1]
point_E = [V_D, 0]
point_A = [V_max, n_max]
point_H = [V_min, -1]
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
plt.plot(V_flaps, n_flaps_list, "black")
plt.plot([V_max_flaps, V_intersect], [2, 2], "black", marker="o")
plt.plot(x_values_AD, y_values_AD, "black", marker="o")
plt.plot(x_values_FE, y_values_FE, 'black', marker="o")
plt.plot(x_values_DE, y_values_DE, "black", marker="o")
plt.plot([point_H[0], point_F[0]], [point_H[1], point_F[1]], "black", marker="o")
plt.plot([0, V_D], [1, 1], 'black', linestyle="--")
plt.plot([V_max, point_A[0]], [0, point_A[1]], 'black', linestyle="--")
plt.plot([V_min, V_min], [0, 1], "black", linestyle="--")
plt.plot([point_F[0], V_cruise], [point_F[1], 0], "black", linestyle="--")
plt.plot([0,0],[0,0],"black",marker="o")
plt.text(0-0.015, 0+0.25, "")
plt.text(point_A[0]-2, point_A[1]+0.2, "A")
plt.text(point_D[0]-2, point_D[1]+0.25, "D")
plt.text(point_E[0]+2, point_E[1]+0.25, "E")
plt.text(point_F[0]-1, point_F[1]-0.25, "F")
plt.text(point_H[0]-2, point_H[1]-0.25, "H")
plt.text(V_min-2,0-0.2, r'$V_{S}$')
plt.text(V_max-2,0-0.2, r'$V_{A}$')
plt.text(V_cruise-2,0+0.1, r'$V_{C}$')
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

#Gust Loads Diagram
CL_alpha = 5.03         #input
chord = b / A                #input
def u_ref(h):
    if h <= 4572:
        u_ref = 17.07 - h * (17.07 - 13.41) / 4572
    elif 4572 < h <= 18288:
        u_ref = 13.41 - (h - 4572) * (13.41 - 6.36) / (18288 - 4572)
    return u_ref

def V_B_calc(h):
    if (h == 0 or h == (5000 * ft_m)):
        T, p, rho, a = ISA_calculator(h, dt_takeoff)
    else:
        T, p, rho, a = ISA_calculator(h, dt_cruise)
    mu = (2 * W_S_design) / (rho * chord * CL_alpha * g)
    Kg = .88 * mu / (5.3 + mu)
    V_B = V_min * np.sqrt(1 + (Kg * rho_0 * u_ref(h) * V_cruise * CL_alpha) / (2 * W_S_design))
    return V_B

def delta_n(h, V):
    if (h == 0 or h == (5000 * ft_m)):
        T, p, rho, a = ISA_calculator(h, dt_takeoff)
    else:
        T, p, rho, a = ISA_calculator(h, dt_cruise)
    delta_n = (rho * V * CL_alpha * u_ref(h)) / (2 * W_S_design)
    return 1 + delta_n, 1 - delta_n
V_B_SL = V_B_calc(0)
delta_cruise_up, delta_cruise_down = delta_n(0, V_cruise)
delta_B_SL_up, delta_B_SL_down = delta_n(0, V_B_SL)
delta_D_up, delta_D_down = delta_n(h_cruise, V_D)
plt.figure()
plt.plot([0, V_B_SL], [1, delta_B_SL_up], "black", linestyle="--", marker="o")
plt.plot([0, V_B_SL], [1, delta_B_SL_down], "black", linestyle="--", marker="o")
plt.plot([0, V_cruise], [1, delta_cruise_up], "black", linestyle="--", marker="o")
plt.plot([0, V_cruise], [1, delta_cruise_down], "black", linestyle="--", marker="o")
plt.plot([0, V_D], [1, delta_D_up], "black", linestyle="--", marker="o")
plt.plot([0, V_D], [1, delta_D_down], "black", linestyle="--", marker="o")
plt.show()


