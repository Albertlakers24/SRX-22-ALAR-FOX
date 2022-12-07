import numpy as np
#Constants
lbs_to_kg = 0.45359237

#Regular efficiencies
wire_eff = 0.99
inverter_eff = 0.995
motor_eff = 0.95
prop_eff = 0.75

#FC efficiencies
FC_eff = 0.6

#Generator efficiencies
turb_eff = 0.47
generator_eff = 0.97

#Total efficiencies
def total_eff(num_inverters, what_eff):
    FC_total_eff = FC_eff * wire_eff * inverter_eff ** num_inverters * motor_eff * prop_eff
    generator_total_eff = turb_eff * generator_eff * wire_eff * inverter_eff ** num_inverters * motor_eff * prop_eff
    if what_eff == 1:
        return FC_total_eff
    if what_eff == 2:
        return generator_total_eff

E_500_nmi = 10354 #MJ
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
max_power_needed = 4032 #kW
FC_mass = max_power_needed / FC_mass_ratio
print(FC_mass)
#Proppeler specs
prop_mass_ratio = 5.22 #kW/kg
prop_mass = max_power_needed / prop_mass_ratio

#1000 nmi mission
num_generators = 2
mass_gen_turb = num_generators * (280 * lbs_to_kg + 335)
E_1000 = 29186

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
perfect_ratio(1)
# for i in np.arange(0, 1.01, 0.01):
#     perfect_ratio(i)

