import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *

T_cruise, p_cruise, rho_cruise = ISA_calculator(h_cruise, 0)
constant_q = 1/2 * CL_max_cruise * rho_cruise
mass = m_mto
S = m_mto * g / W_S_design
V = []
n_list = []
for i in np.arange(0, V_cruise, 0.1):
    V.append(i)
    n = constant_q * i**2 / W_S_design
    n_list.append(n)
print(n_list)
plt.plot(V, n_list)
plt.show()