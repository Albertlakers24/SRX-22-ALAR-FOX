#Constants
g = 9.80665                         #gravity
e_lh2 = 120                         #MJ/kg Specific Energy Liquid Hydrogen
kts_m_s = 0.514444                  # knots to m/s
ft_m = 0.3048                       #ft to m
FL_ft = 100                         #FL to ft
V_cruise = 275 * kts_m_s            #Cruise Velocity Requirement
h_cruise = 280 * FL_ft * ft_m       #Cruise height Requirement
PAX = 48                            #Passenger Requirement
T_0 = 288.15                        #ISA Temperature
p_0 = 101325                        #ISA Pressure
rho_0 = 1.225                       #ISA Density
lapse_rate = -6.5/1000              #Troposphere lapse rate
Molar_mass_air = 8.31432            #N m kmol⁻¹ K⁻¹
specific_gas_constant = 287.052     #J·kg⁻¹·K⁻¹
universal_gas_constant = 8.31432    #N m kmol⁻¹ K⁻¹





#General Functions
def ISA_calculator(h):
    T = T_0 + lapse_rate * h
    p = p_0 * ((T_0 / T) ** ((g * Molar_mass_air) / (universal_gas_constant * lapse_rate)))
    rho = p / (specific_gas_constant * T)
    return T, p, rho
