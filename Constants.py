#Constants
g = 9.80665                         #gravity
e_lh2 = 120*10**6                   #J/kg Specific Energy Liquid Hydrogen
kts_m_s = 0.514444                  #knots to m/s
nmi_m = 1852                        #nmi to m
min_s = 60                          #minutes to seconds
ft_m = 0.3048                       #ft to m
km_m = 1000                         #km to m
FL_ft = 100                         #FL to ft
lbs_kg = 0.453592                   #lbs to kg
V_cruise = 275 * kts_m_s            #Cruise Velocity Requirement
h_cruise = 280 * FL_ft * ft_m       #Cruise height Requirement
PAX = 48                            #Passenger Requirement
T_0 = 288.15                        #ISA Temperature
p_0 = 101325                        #ISA Pressure
rho_0 = 1.225                       #ISA Density
lapse_rate = -6.5/km_m              #Troposphere lapse rate
Molar_mass_air = 8.31432            #N m kmol⁻¹ K⁻¹
specific_gas_constant = 287.052     #J·kg⁻¹·K⁻¹
universal_gas_constant = 8.31432    #N m kmol⁻¹ K⁻¹
E = 30 * 60                         #Loiter endurance in seconds
s_takeoff = 4500 * ft_m             #Takeoff Distance
s_landing = 4500 * ft_m             #Landing Distance
CL_max_take0ff = 2.1                # -
CL_max_cruise = 1.9                 # -
CL_max_landing = 2.6                # -
WPAX = 200*lbs_kg*PAX*g             #Fat American Weight
WPAXBAGGAGE = 40*lbs_kg*PAX*g       #Baggage Weight for Passenger
R_norm = 1000 * nmi_m               #Design Range
R_div = 100 * nmi_m                 #Divergence Range
Loiter = 30 * min_s                 #Loiter Endurance
eta_prop = 0.85                     #Propeller efficiency
eta_EM = 0.95                       #Electric motor efficiency
eta_wire = 0.97                     #Wire efficiency
eta_inverter = 0.995                #Inverter efficiency
eta_fuelcell = 0.60                 #Fuel cell efficiency

#General Functions
def ISA_calculator(h):
    T = T_0 + lapse_rate * h
    p = p_0 * ((T_0 / T) ** ((g * Molar_mass_air) / (universal_gas_constant * lapse_rate)))
    rho = p / (specific_gas_constant * T)
    return T, p, rho
