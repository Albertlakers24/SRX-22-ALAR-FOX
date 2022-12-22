import numpy as np
import matplotlib.pyplot as plt
V_cruise = 275 #KTAS
kts_to_ms = 0.51444444444444
ft_to_m = 0.3048
V_cruise_ms = V_cruise * kts_to_ms
M_tip = 0.95

E_500_nmi = 11835 #MJ
E_1000 = 34943 #MJ
max_power_needed_all_eff = 3835 #kW
#Regular efficiencies
wire_eff = 0.97
inverter_eff = 0.995
motor_eff = 0.95
prop_eff = 0.85

#FC efficiencies
FC_eff = 0.6

max_power_needed_no_eff = max_power_needed_all_eff * FC_eff * wire_eff * inverter_eff**2 * motor_eff * prop_eff
max_power_needed = max_power_needed_no_eff / motor_eff / prop_eff
print(max_power_needed)
power_cruise = 3148 #2758 #kW
max_power_needed = 4669#4946#2474

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
V_max_cruise = M_tip * a_cruise
V_takeoff = 65 #m/s
V_cruise = 275 * kts_to_ms
def tip_speed(V_real, V_max):
    V_tip = np.sqrt(V_max**2 - V_real**2)
    return V_tip
V_tip_takeoff = tip_speed(V_takeoff, V_max_takeoff)
V_tip_cruise = tip_speed(V_cruise, V_max_cruise)

def diameter(prop_blades, power, mass):
    if prop_blades == 2:
        D = 0.56 * (power * mass) ** (1/4)
    elif prop_blades == 3:
        D = 0.52 * (power * mass) ** (1/4)
    elif prop_blades == 4:
        D = 0.49 * (power * mass) ** (1/4)
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
ratio_big_small = 0.8 #power_cruise / max_power_needed
power_big = max_power_needed * ratio_big_small / number_big
mass_big = power_big / power
big_prop_diameter = diameter(3, power, mass_big)
rpm_big_takeoff = rpm(big_prop_diameter, V_tip_takeoff)
rpm_big_cruise = rpm(big_prop_diameter, V_tip_cruise)
print(big_prop_diameter, rpm_big_takeoff, rpm_big_cruise)

def number_props(num_blades, ratio):
    diameters_small = []
    number_props = []
    rpms_cruise = []
    rpms_takeoff = []
    for j in range(2, 15, 2):
        number_props.append(j)
        power_small = max_power_needed * (1 - ratio) / j
        mass_small_est_val = power_small / power
        new_diameter = diameter(num_blades, power, mass_small_est_val)
        diameters_small.append(new_diameter)
        rpms_cruise.append(rpm(new_diameter, V_tip_cruise))
        rpms_takeoff.append(rpm(new_diameter, V_tip_takeoff))
    return diameters_small, number_props, rpms_cruise, rpms_takeoff


# power_big = max_power_needed * ratio_big_small / number_big
# mass_big_est = power_big / power
# mass_small = power1 / power
# mass_big = power2 / power
# diameters = [diameter(2, power, mass_small), diameter(3, power, mass_small), diameter(2, power, mass_big), diameter(3, power, mass_big)]
#
# rpms = []
# for i in np.arange(0, len(diameters)):
#     rpms.append(rpm(diameters[i], V_tip_cruise))
#     rpms.append(rpm(diameters[i], V_tip_takeoff))

# print(f"With a power of {power1} kW, the two blade propeller has a diameter of {diameters[0]} m and the three blade propeller of {diameters[1]}")
# print(f"With a power of {power2} kW, the two blade propeller has a diameter of {diameters[2]} m and the three blade propeller of {diameters[3]}")

diameters_small_2, number_props_2, rpms_cruise_2, rpms_takeoff_2 = number_props(2, ratio_big_small)
diameters_small_3, number_props_3, rpms_cruise_3, rpms_takeoff_3 = number_props(3, ratio_big_small)
diameters_small_4, number_props_4, rpms_cruise_4, rpms_takeoff_4 = number_props(4, ratio_big_small)
print(diameters_small_3)
# print(rpms_cruise, rpms_takeoff)
# diameters_small_split = np.array_split(diameters_small, 10)
# number_props_split = np.array_split(number_props, 10)
diameters_accumulated_2 = np.array(diameters_small_2) * np.array(number_props_2)
diameters_accumulated_3 = np.array(diameters_small_3) * np.array(number_props_3)
diameters_accumulated_4 = np.array(diameters_small_4) * np.array(number_props_4)
total_diameter_3 = []
total_diameter_2 = []
total_diameter_4 = []
fuselage_width = 2.8
for i in np.arange(0, 7, 1):
    total_diameter_3.append(diameters_accumulated_3[i] + big_prop_diameter + 0.04 * diameters_small_3[i] * 2 * (number_props_3[i]) + fuselage_width)
    total_diameter_2.append(diameters_accumulated_2[i] + big_prop_diameter + 0.04 * diameters_small_2[i] * 2 * (number_props_2[i]) + fuselage_width)
    total_diameter_4.append(diameters_accumulated_4[i] + big_prop_diameter + 0.04 * diameters_small_4[i] * 2 * (number_props_4[i]) + fuselage_width)
print(total_diameter_4)
print(rpms_takeoff_3)


plt.figure("Number of props vs diameter per prop, ratio 0.6")
plt.title("Number of props vs diameter per propeller")
plt.xlabel("Number of propellers")
plt.ylabel("Propeller diameter (m)")
plt.plot(number_props_2, diameters_small_2, "blue", label = "2 blades")
plt.plot(number_props_3, diameters_small_3, "red", label = "3 blades")
plt.plot(number_props_4, diameters_small_4, "green", label = "4+ blades")
plt.legend()
# plt.savefig('ratio 0.5.png')
plt.show()
plt.figure("Number of props vs accumulated diameters, ratio 0.6")
plt.title("Number of props vs accumulated diameters")
plt.xlabel("Number of propellers")
plt.ylabel("Propeller diameters accumulated (m)")
plt.plot(number_props_2, diameters_accumulated_2, "blue", label = "2 blades")
plt.plot(number_props_3, diameters_accumulated_3, "red", label = "3 blades")
plt.plot(number_props_4, diameters_accumulated_4, "green", label = "4+ blades")
plt.legend()
# plt.savefig("ratio 0.5 acc.png")
plt.show()
plt.figure("props_vs_total")
plt.title("Total distance taken on wing vs number of propellers")
plt.xlabel("Number of propellers")
plt.ylabel("Total distance taken by propellers")
plt.plot(number_props_2, total_diameter_2, "blue", label = "2 blades")
plt.plot(number_props_3, total_diameter_3, "red", label = "3 blades")
plt.plot(number_props_4, total_diameter_4, "green", label = "4+ blades")
plt.legend()
plt.show()

# del rpms_cruise_3[6]
# del rpms_cruise_3[5]
num_props_smaller_3 = number_props_3.copy()
# del num_props_smaller_3[6]
# del num_props_smaller_3[5]
# del rpms_cruise_2[6]
# del rpms_cruise_2[5]
num_props_smaller_2 = number_props_2.copy()
# del num_props_smaller_2[6]
# del num_props_smaller_2[5]
# del rpms_cruise_4[6]
# del rpms_cruise_4[5]
num_props_smaller_4 = number_props_4.copy()
# del num_props_smaller_4[6]
# del num_props_smaller_4[5]
plt.figure("RPM vs number of propellers")
plt.plot(num_props_smaller_2, rpms_cruise_2, "blue", label = "2 blades")
plt.plot(num_props_smaller_3, rpms_cruise_3, "red", label = "3 blades")
plt.plot(num_props_smaller_4, rpms_cruise_4, "green", label = "4+ blades")
plt.legend()
plt.title("RPM vs number of propellers")
plt.xlabel("Number of propellers")
plt.ylabel("RPM of small propellers during take-off")
plt.show()
