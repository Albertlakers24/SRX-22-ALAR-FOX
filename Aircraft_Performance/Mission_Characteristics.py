import numpy as np
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto, m_f
from Constants import *
from Initial_Aircraft_Sizing.Wing_planform import Sw
from Initial_Aircraft_Sizing.Propeller_Selection import *

s_0 = 0
s_1 = 5000 * ft_m
s_2 = 15000 * ft_m
s_3 = 28000 * ft_m
s_4 = 10000 * ft_m
s_5 = 0
ROC1 = 4
ROC2 = 3
ROC3 = 2
ROD1 = -1500 * ft_m / 60
ROD2 = -1110 * ft_m / 60

V_IAS1 = 140 * kts_m_s
V_IAS2 = 176 * kts_m_s
V_IAS3 = 176 * kts_m_s
V_IAS4 = 176 * kts_m_s
V_IAS5 = 140 * kts_m_s
V_IAS_list = [V_IAS1, V_IAS2, V_IAS3]
s_list = [s_0, s_1, s_2, s_3, s_4, s_5]
V_TAS_list = []
def V_TAS_calc(V_IAS, h, h_0):
    V_TAS = V_IAS * (1 + (h + h_0)/ 2 / 1000 * 0.02)
    V_TAS_list.append(V_TAS)
    return V_TAS

V_TAS1_avg = V_TAS_calc(V_IAS1, s_1, s_0)
V_TAS2_avg = V_TAS_calc(V_IAS2, s_2, s_1)
V_TAS3_avg = V_TAS_calc(V_IAS3, s_3, s_2)
V_TAS4_avg = V_TAS_calc(V_IAS4, s_3, s_4)
V_TAS5_avg = V_TAS_calc(V_IAS5, s_4, s_5)
# print(V_TAS_list)
ROC_list = [ROC1, ROC2, ROC3, ROD1, ROD2]

V_x = []
for i in np.arange(0, len(V_TAS_list), 1):
    V_x_part = np.sqrt(V_TAS_list[i]**2 - ROC_list[i]**2)
    V_x.append(V_x_part)

def t_calc():
    t_list = []
    for i in np.arange(0, 3):
        t = (s_list[i+1] - s_list[i]) / ROC_list[i]
        t_list.append(t)
    for i in np.arange(3,5):
        t = (s_list[i] - s_list[i+1]) / ROC_list[i] * -1
        t_list.append(t)
    return t_list
t_list = t_calc()
# print("times", t_list)

total_distance = sum(np.array(V_x) * np.array(t_list))
# print(total_distance / 1852)
full_time_climb = t_list[0] + t_list[1] + t_list[2]
full_time_descent = t_list[3] + t_list[4]
# print(f"Full time of climb: {full_time_climb} s")
# print(f"Full time of descent: {full_time_descent} s")
cruise_distance = (1000 * 1852 - total_distance)                             #in m
cruise_time = cruise_distance / V_cruise
# print(f"Cruise time: {cruise_time} s")
full_time = full_time_climb + cruise_time + full_time_descent
# print(f"Full time of 1000 nmi mission: {full_time} s")

fuel_flow = m_f / full_time                                                     #Fuel Flow per second
# print(fuel_flow)

Weight = m_mto
C_L = 0.72
Lift = 1/2 * rho_cruise * Sw * C_L * V_cruise**2
CL_need = (Weight * 9.80665) / (1/2 * rho_cruise * Sw * V_cruise**2)
S_opt = (Weight * 9.80665) / (1/2 * rho_cruise * C_L * V_cruise**2)

# print(Lift - Weight * 9.80665)
# print(CL_need)
# print(S_opt)
# print(np.sqrt(S_opt * 12))

# THRUST AND POWER CALCULATIONS FOR PROPELLERS
S_big_prop = np.pi / 4 * big_prop_diameter ** 2
CL_here = 0.72
CD_here = CL_here / max_CL_CD


# V_takeoff can be used as V_0
def V_e(ROC, V, W):
    flight_angle = np.arcsin(ROC / V)
    D = rho_5000 * CD_here * V ** 2 * Sw / 2
    T = W * np.sin(flight_angle) + D
    Ve = np.sqrt(T / (0.5 * rho_5000 * S_big_prop) + V ** 2)
    return Ve


def mass_flow(V, h_start, h_end, dt, CL, ROC):
    h = (h_start + h_end) / 2
    T, p, rho, a = ISA_calculator(h, dt)
    CD = CL / max_CL_CD
    D = 0.5 * rho * V ** 2 * Sw * CD
    Power_req = D * V
    Power_available = ROC * m_mto + Power_req
    mass_flow = Power_available / (H2_power_density * 10**6)
    return mass_flow

climb1 = mass_flow(V_TAS1_avg, s_0, s_1, 0, CL_here, ROC1)
climb2 = mass_flow(V_TAS2_avg, s_1, s_2, 0, CL_here, ROC2)
climb3 = mass_flow(V_TAS3_avg, s_2, s_3, 0, CL_here, ROC3)
cruise = mass_flow(V_cruise, s_3, s_3, 0, CL_here, 0)
descent1 = mass_flow(V_TAS4_avg, s_3, s_4, 0, CL_here, ROD1)
descent2 = mass_flow(V_TAS5_avg, s_4, s_5, 0, CL_here, ROD2)
fuel_flows = [climb1, climb2, climb3, cruise, descent1, descent2]
total_fuel = 327
total_fuel_flows = sum(fuel_flows)
ff_1 = climb1 / climb1
ff_2 = climb2 / climb1
ff_3 = climb3 / climb1
ff_4 = cruise / climb1
ff_5 = descent1 / climb1
ff_6 = descent2 / climb1
ffs = [ff_1, ff_2, ff_3, ff_4, ff_5, ff_6]

t_list_all = [t_list[0], t_list[1], t_list[2], cruise_time, t_list[3], t_list[4]]

first_ff = total_fuel / sum(np.array(t_list_all) * np.array(ffs))
second_ff = ff_2 * first_ff
third_ff = ff_3 * first_ff
fourth_ff = ff_4 * first_ff
fifth_ff = ff_5 * first_ff
sixth_ff = ff_6 * first_ff
all_ffs = [first_ff, second_ff, third_ff, fourth_ff, fifth_ff, sixth_ff]
