import numpy as np
import matplotlib.pyplot as plt
V_cruise = 275 #KTAS
kts_to_ms = 0.51444444444444
ft_to_m = 0.3048
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
    a = np.sqrt(gamma * T * specific_gas_constant)
    return T, rho, a
T_sea, rho_sea, a_sea = ISA_calculator(0 * ft_to_m)
T_cruise, rho_cruise, a_cruise = ISA_calculator(28000 * ft_to_m)
M_cruise = V_cruise_ms / a_cruise
rps = 2500 / 60
# D = a / (np.pi * rps) * np.sqrt(M_tip**2 - M_cruise**2)
# print(rho, a)

V_max_takeoff = 213 #m/s
V_max_cruise = 0.8 * a_cruise
V_takeoff = 65 #m/s
V_cruise = 275 * kts_to_ms
def tip_speed(V_real, V_max):
    V_tip = np.sqrt(V_max**2 - V_real**2)
    return V_tip
V_tip_takeoff = tip_speed(V_takeoff, V_max_takeoff)
V_tip_cruise = tip_speed(V_cruise, V_max_cruise)

def diameter(prop_blades, power, mass):
    if prop_blades == 2:
        D = 0.6 * (power * mass) ** (1/4)
    elif prop_blades == 3:
        D = 0.5 * (power * mass) ** (1/4)
    return D

power1 = 190 #kW
power2 = 2660 #kW


def rpm(D, V_tip):
    rpm = V_tip / (np.pi * D) * 60
    return rpm
    
#Method 2
power = 15 #kW / kg
number_big = 2
number_small = 10
ratio_big_small = 0.5

def number_props(ratio):
    diameters_small = []
    number_props = []
    # for i in np.arange(0, 1.1, 0.1):
    #     power_big = max_power_needed * i / number_big
    #     mass_big = power_big / power
    #     diameter_big = diameter(3, power, mass_big)
    for j in range(2, 15, 2):
        number_props.append(j)
        power_small = max_power_needed * (1 - ratio) / j
        mass_small_est_val = power_small / power
        diameters_small.append(diameter(3, power, mass_small_est_val))
    return diameters_small, number_props


power_big = max_power_needed * ratio_big_small / number_big
mass_big_est = power_big / power
# mass_small_est = []
# total_mass_small = []
# def number_props(ratio):
#     diameters_small = []
#     number_props = []
#     for i in np.arange(2, 15, 2):
#         number_props.append(i)
#         power_small = max_power_needed * (1 - ratio) / i
#         mass_small_est_val = power_small / power
#         mass_small_est.append(mass_small_est_val)
#         diameters_small.append(diameter(3, power, mass_small_est_val))
#     return diameters_small, number_props
# print("Number of props", number_props(ratio_big_small))
mass_small = power1 / power
mass_big = power2 / power
diameters = [diameter(2, power, mass_small), diameter(3, power, mass_small), diameter(2, power, mass_big), diameter(3, power, mass_big)]
rpms = []
for i in np.arange(0, len(diameters)):
    rpms.append(rpm(diameters[i], V_tip_cruise))
    rpms.append(rpm(diameters[i], V_tip_takeoff))

print(f"With a power of {power1} kW, the two blade propeller has a diameter of {diameters[0]} m and the three blade propeller of {diameters[1]}")
print(f"With a power of {power2} kW, the two blade propeller has a diameter of {diameters[2]} m and the three blade propeller of {diameters[3]}")
print(rpms)
print(diameter(3, power, mass_big))
diameters_small, number_props = number_props(ratio_big_small)
# diameters_small_split = np.array_split(diameters_small, 10)
# number_props_split = np.array_split(number_props, 10)
diameters_accumulated = np.array(diameters_small) * np.array(number_props)
plt.figure("Number of props vs diameter per prop, ratio 0.5")
plt.plot(number_props, diameters_small)
plt.savefig('ratio 0.5.png')
plt.show()
plt.figure("Number of props vs accumulated diameters, ratio 0.5")
plt.plot(number_props, diameters_accumulated)
plt.savefig("ratio 0.5 acc.png")
plt.show()
