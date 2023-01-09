import numpy as np
#Inputs
E_500_nmi = 11835 #MJ
E_1000 = 34943 #MJ
max_power_needed = 3760 #kW
ratio = 1 #Ratio of H2 / fuel

#Constants
lbs_to_kg = 0.45359237

#Regular efficiencies
wire_eff = 0.97
inverter_eff = 0.995
motor_eff = 0.95
prop_eff = 0.85

#FC efficiencies
FC_eff = 0.6

#Generator efficiencies
turb_eff = 0.39
generator_eff = 0.97

#Total efficiencies
def total_eff(num_inverters, what_eff):
    FC_total_eff = FC_eff * wire_eff * inverter_eff ** num_inverters * motor_eff * prop_eff
    generator_total_eff = turb_eff * generator_eff * wire_eff * inverter_eff ** num_inverters * motor_eff * prop_eff
    if what_eff == 1:
        return FC_total_eff
    if what_eff == 2:
        return generator_total_eff

E_500_req = E_500_nmi / total_eff(2, 1)
H2_E_density = 120 #MJ/kg
mass_H2_500 = E_500_req / H2_E_density
H2_density = 66.31 #kg/m3
volume_H2_500 = mass_H2_500 / H2_density
ratio_tank = 67 / 150
mass_tanks = mass_H2_500 * ratio_tank
mass_H2_total_500 = mass_tanks + mass_H2_500
print(f"For the 500 nmi, mission, we need {mass_H2_total_500} kg and {volume_H2_500} m3 of H2")

#FC mass
FC_mass_ratio = 3 #kW/kg
FC_mass = max_power_needed / FC_mass_ratio
print(FC_mass)
#Proppeler specs
prop_mass_ratio = 5.22 #kW/kg
prop_mass = max_power_needed / prop_mass_ratio

#Inverter mass
power_density_inverter = 30 #kW/kg
mass_inverter = max_power_needed / power_density_inverter
print(mass_inverter)

#1000 nmi mission



#Specs jet fuel
jet_E_density = 42.8 #MJ/kg
jet_density = 820 #kg/m3
volume_H2_500 = mass_H2_500 / H2_density
ratio_tank = 67 / 150
mass_tanks = mass_H2_500 * ratio_tank

#Give ratio of H2 / fuel
def perfect_ratio(ratio):
    E_H2 = ratio * E_1000 / total_eff(2, 1)
    E_Fuel = (1 - ratio) * E_1000 / total_eff(2, 2)
    mass_H2 = E_H2 / H2_E_density
    mass_fuel = E_Fuel / jet_E_density
    volume_H2 = mass_H2 / H2_density
    volume_fuel = mass_fuel / jet_density
    mass_H2_tank = mass_H2 * ratio_tank
    mass_H2_total = mass_H2 + mass_H2_tank
    print(f"Mass fuel is {mass_fuel}, volume is {volume_fuel}",
         f"Mass H2 + tanks is {mass_H2_total}, volume is {volume_H2}")
    total_fuel_mass = mass_fuel + mass_H2_total
    return total_fuel_mass
# perfect_ratio(1)
# for i in np.arange(0, 1.01, 0.01):
#     perfect_ratio(i)
if ratio == 1:
    num_generators = 0
elif ratio == 0.75:
    num_generators = 1
elif ratio == 0.5:
    num_generators = 2
elif ratio == 0.75:
    num_generators = 3
elif ratio == 0:
    num_generators = 4
    FC_mass = 0
mass_gen_turb = num_generators * (280 * lbs_to_kg + 335)
total_weight = FC_mass + prop_mass + mass_gen_turb + mass_inverter + perfect_ratio(ratio)
print(total_weight)
print(total_eff(2, 1))
# masses = []
# for ratio in np.arange(1, -0.25, -0.25):
#     if ratio == 1:
#         num_generators = 0
#     elif ratio == 0.75:
#         num_generators = 1
#     elif ratio == 0.5:
#         num_generators = 2
#     elif ratio == 0.75:
#         num_generators = 3
#     elif ratio == 0:
#         num_generators = 4
#         FC_mass = 0
#     mass_gen_turb = num_generators * (280 * lbs_to_kg + 335)
#     total_weight = FC_mass + prop_mass + mass_gen_turb + perfect_ratio(ratio) + mass_inverter
#     masses.append(total_weight)
# print(masses)
# print(total_weight)


# total_weight_1 = FC_mass + prop_mass + mass_gen_turb + perfect_ratio(1) + mass_inverter
# total_weight_75 = FC_mass + prop_mass + mass_gen_turb + perfect_ratio(0.75) + mass_inverter
# total_weight_50 = FC_mass + prop_mass + mass_gen_turb + perfect_ratio(0.5) + mass_inverter
# total_weight_25 = FC_mass + prop_mass + mass_gen_turb + perfect_ratio(0.25) + mass_inverter
# total_weight_0 = FC_mass + prop_mass + mass_gen_turb + perfect_ratio(0) + mass_inverter
# print(total_weight_1, total_weight_75, total_weight_50, total_weight_25, total_weight_0)