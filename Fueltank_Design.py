import numpy as np

operating_design_stress = 172 #MPa
ultimate_design_stress = 234 #MPa
allowable_stress = 30000 #psi

def stress(min_max_stress_ratio, R2):
    sigma_a1 = operating_design_stress / (1 - (operating_design_stress * 0.5 * (1 + min_max_stress_ratio)) / ultimate_design_stress)
    stress = sigma_a1 / (1 + (sigma_a1 * 0.5 * (1 + min_max_stress_ratio)) / ultimate_design_stress)
    return stress

P_vent = 0.300000 #MPa
r_shell = 2.5 / 2 #m
safety_factor = 0.8

#Constants
p_0 = 0.101325 #MPa
T_0 = 288.15 #K
lapse_rate = -0.0065 #K/m
g_0 = 9.80665
Molar_mass_air = 0.0289644 #kg/mol
universal_gas_constant = 8.31432 #N m kmol⁻¹ K⁻¹

def t_shell_central(h):
    T = T_0 + lapse_rate * h
    P_out = p_0 * ((T_0 / T) ** ((g_0 * Molar_mass_air) / (universal_gas_constant * lapse_rate)))
    thickness_shell_central = ((P_vent - P_out) * r_shell) / (stress(0.43, 0.43) * safety_factor)
    return thickness_shell_central

print(t_shell_central(28000 * 0.3048))

