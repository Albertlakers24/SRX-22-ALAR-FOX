import numpy as np
V_cruise = 275 #KTAS
kts_to_ms = 0.51444444444444
V_cruise_ms = V_cruise * kts_to_ms
M_tip = 0.8

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
T, rho = ISA_calculator(0 * 0.3048)
a = np.sqrt(gamma * T * specific_gas_constant)
M_cruise = V_cruise_ms / a
rps = 2500 / 60
D = a / (np.pi * rps) * np.sqrt(M_tip**2 - M_cruise**2)
print(rho, a)

