import numpy as np
from Class_I_Weight_Estimation import Cd0, A, e, MTOW

#Constants
g_0 = 9.80665
Molar_mass_air = 0.0289644 #kg/mol
universal_gas_constant = 8.31432 #N m kmol⁻¹ K⁻¹
specific_gas_constant = 287.052 #J·kg⁻¹·K⁻¹

b = 25 #wingspan in m
S_w = b**2 / A
kts_ms = 0.514444444
V_cruise = 275 * kts_ms #KTAS
T_0 = 288.15
p_0 = 101325 #Pa
lapse_rate = -0.0065


def ISA_calculator(h):
    T = T_0 + lapse_rate * h
    p = p_0 * ((T_0 / T) ** ((g_0 * Molar_mass_air) / (universal_gas_constant * lapse_rate)))
    rho = p / (specific_gas_constant * T)
    return T, p, rho

T, p, rho = ISA_calculator(28000 * 0.3048)

def C_L_calc(W, rho, V):
    C_L_calc = 2 * W / (rho * V**2 * S_w)
    return C_L_calc

def C_D_calc(AR, W, rho, V, C_D_0):
    C_D = C_D_0 + C_L_calc(W, rho, V)**2 / (np.pi * e * AR)
    return C_D, rho

print(C_D_calc(A, MTOW, rho, V_cruise, Cd0))
print(C_L_calc(MTOW, rho, V_cruise))
