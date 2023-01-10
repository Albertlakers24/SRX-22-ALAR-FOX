import numpy as np
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto, m_f
from Constants import *
from Initial_Aircraft_Sizing.Wing_planform import Sw

s_0 = 0
s_1 = 5000 * ft_m
s_2 = 15000 * ft_m
s_3 = 28000 * ft_m
ROC1 = 4
ROC2 = 3
ROC3 = 2


V_IAS1 = 140 * kts_m_s
V_IAS2 = 176 * kts_m_s
V_IAS3 = 176 * kts_m_s
V_IAS_list = [V_IAS1, V_IAS2, V_IAS3]
s_list = [s_0, s_1, s_2, s_3]
V_TAS_list = []
def V_TAS_calc(V_IAS, h, h_0):
    V_TAS = V_IAS * (1 + (h + h_0)/ 2 / 1000 * 0.02)
    V_TAS_list.append(V_TAS)
    return V_TAS

V_TAS1_avg = V_TAS_calc(V_IAS1, s_1, s_0)
V_TAS2_avg = V_TAS_calc(V_IAS2, s_2, s_1)
V_TAS3_avg = V_TAS_calc(V_IAS3, s_3, s_2)
print(V_TAS_list)
ROC_list = [ROC1, ROC2, ROC3]

V_x = []
for i in np.arange(0, len(V_IAS_list), 1):
    V_x_part = np.sqrt(V_TAS_list[i]**2 - ROC_list[i]**2)
    V_x.append(V_x_part)

def t_calc():
    t_list = []
    for i in np.arange(0, 3):
        t = (s_list[i+1] - s_list[i]) / ROC_list[i]
        t_list.append(t)
    return t_list
t_list = t_calc()
print(t_list)

total_distance = sum(np.array(V_x) * np.array(t_list))
print(V_x)
print(total_distance / 1852)
full_time_climb = sum(t_list)
print(f"Full time of climb and descent: {full_time_climb} s")
remaining_distance = (1000 * 1852 - total_distance)                             #in m
full_time = full_time_climb + remaining_distance / V_cruise
print(f"Full time of 1000 nmi mission: {full_time} s")

fuel_flow = m_f / full_time                                                     #Fuel Flow per second
print(fuel_flow)

Weight = m_mto
C_L = 0.72
Lift = 1/2 * rho_cruise * Sw * C_L * V_cruise**2
CL_need = (Weight * 9.80665) / (1/2 * rho_cruise * Sw * V_cruise**2)
S_opt = (Weight * 9.80665) / (1/2 * rho_cruise * C_L * V_cruise**2)

# print(Lift - Weight * 9.80665)
# print(CL_need)
# print(S_opt)
# print(np.sqrt(S_opt * 12))
