import numpy as np
# import class i weight estimation

S_w = 5
Surfaces_ratio = 6
C_fe = 0.5
e = 0.3
rho_sea = 1.225 #kg/m3
rho_cruise = 0.02047 * rho_sea
V_cruise = 275 #KTAS
h_cruise = 28000 #ft
ft_to_m = 0.3048

C_D_0 = Surfaces_ratio * C_fe

def C_L_calc(W, rho, V):
    C_L_calc = 2 * W / (rho * V**2 * S_w)
    return C_L_calc

def C_D_calc(AR, W, rho, V):
    C_D = C_D_0 + C_L_calc(W, rho, V)**2 / (np.pi * e * AR)
    return C_D

print(C_D_calc(5, 23405, 1.0, 300))


