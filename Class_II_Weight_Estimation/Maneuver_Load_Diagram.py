import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *

T_cruise, p_cruise, rho_cruise = ISA_calculator(h_cruise, 0)
constant_q = 1/2 * CL_max_cruise * rho_cruise
mass = m_mto
S = m_mto * g / W_S_design
V_A = []
n_A_list = []
n_max = max(2.1 * (24000 / (m_mto + 10000)), 2.5)
n_min = -1
V_max = np.sqrt((2 * n_max * W_S_design) / (rho_cruise * CL_max_cruise))
for i in np.arange(0, V_max, 0.1):
    V_A.append(i)
    n = constant_q * i**2 / W_S_design
    n_A_list.append(n)

V_H = []
n_H_list = []
V_min = np.sqrt((2 * 1 * W_S_design) / (rho_cruise * CL_max_cruise))
for j in np.arange(0, V_min, 0.1):
    V_H.append(j)
    n = -constant_q * j ** 2 / W_S_design
    n_H_list.append(n)

n_F_list = []
V_F = []
for i in np.arange(V_min, V_cruise, 0.1):
    n_F_list.append(-1)
    V_F.append(i)

V_D = V_cruise / 0.8
V_D_list = []
n_D_list = []
for i in np.arange(V_max, V_D, 0.1):
    n_D_list.append(n_max)
    V_D_list.append(i)

plt.figure()
plt.plot(V_A, n_A_list, "blue")
plt.plot(V_H, n_H_list, "blue")
plt.plot(V_F, n_F_list, "blue")
plt.plot(V_D_list, n_D_list, "blue")
plt.xlabel("V")
plt.ylabel("n")

plt.show()