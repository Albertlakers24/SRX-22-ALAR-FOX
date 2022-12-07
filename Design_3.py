import numpy as np
#Constants
lbs_to_kg = 0.45359237

#Regular efficiencies
wire_eff = 0.99
inverter_eff = 0.995
motor_eff = 0.95
prop_eff = 0.75

#FC efficiencies
FC_eff = 0.65

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

print(total_eff(2, 1), total_eff(2, 2))

E_500_nmi = 10354 #MJ
E_500_req = E_500_nmi / total_eff(2, 1)
H2_E_density = 120 #MJ/kg
mass_H2_500 = E_500_req / H2_E_density
H2_density = 66.31 #kg/m3
volume_H2_500 = mass_H2_500 / H2_density
ratio_tank = 67 / 150
mass_tanks = mass_H2_500 * ratio_tank
print(mass_tanks)

#FC mass
FC_mass_ratio = 11 #kW/kg
max_power_needed = 4032 #kW
FC_mass = max_power_needed / FC_mass_ratio

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
print(mass_tanks)

#Give ratio of H2 / fuel
def perfect_ratio(ratio):
    E_H2 = ratio * E_1000 / total_eff(2, 1)
    E_Fuel = (1 - ratio) * E_1000 / total_eff(2, 2)
    mass_H2 = E_H2 / H2_E_density
    mass_fuel = E_Fuel / jet_E_density



