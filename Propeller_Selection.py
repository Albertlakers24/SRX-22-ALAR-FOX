import numpy as np
V_cruise = 275 #KTAS
kts_to_ms = 0.51444444444444
V_cruise_ms = V_cruise * kts_to_ms
M_tip = 0.8

E_500_nmi = 11835 #MJ
E_1000 = 34943 #MJ
max_power_needed = 3760 #kW
power_cruise = 2758 #kW

g_0 = 9.80665
Molar_mass_air = 0.0289644 #kg/mol
universal_gas_constant = 8.31432 #N m kmol⁻¹ K⁻¹
specific_gas_constant = 287.052 #J·kg⁻¹·K⁻¹
gamma = 1.4

T_0 = 288.15
p_0 = 101325 #Pa
lapse_rate = -0.0065

def ISA_calculator(h):
    T = T_0 + lapse_rate * h
    p = p_0 * ((T_0 / T) ** ((g_0 * Molar_mass_air) / (universal_gas_constant * lapse_rate)))
    rho = p / (specific_gas_constant * T)
    return T, rho
T_sea, rho_sea = ISA_calculator(0 * 0.3048)
a_sea = np.sqrt(gamma * T_sea * specific_gas_constant)
# M_cruise = V_cruise_ms / a
rps = 2500 / 60
# D = a / (np.pi * rps) * np.sqrt(M_tip**2 - M_cruise**2)
# print(rho, a)

V_max = 213 #m/s
V_takeoff = 65 #m/s
def tip_speed(V_real):
    V_tip = np.sqrt(V_max**2 - V_real**2)
    return V_tip
print(tip_speed(V_takeoff))

