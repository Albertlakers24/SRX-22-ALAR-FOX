import numpy as np
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import *

cg_start = OEW_cg                                       #OEW center of gravity from nose (m)
cargo_front = 123                                       #Front cargo position from the nose (m)
cargo_aft = 13                                          #Aft cargo position from the nose (m)
l_cabin = 13.0                                          #Cabin length (m)
overhead_surface = 0.1508                               #Surface area per overhead bin (m^2)
overhead_volume = overhead_surface * 2 * l_cabin        #Overhead volume available
req_volume = vol_pax_baggage                            #Cargo volume required
front_volume = max(0, req_volume - overhead_volume)     #Cargo volume in front
front_surface_area = 1.37                               #Cargo surface area per front compartment (m^2)
len_front_storage = front_volume/(2*front_surface_area) #Length of front compartment (m)
mass_front = front_volume/vol_pax_baggage*m_pax_baggage #Mass front storage compartment (kg)
mass_overhead = m_pax_baggage - mass_front              #Mass overhead storage (kg)
x_front = 3.3                                           #Start of circular part (m)
x_cg_front = 3.3 + len_front_storage / 2                #CG location of front storage (m)
#POINT LOADS APPLY AT EVERY SEAT, SO 13/12, / 2 to obtain its cg
x_cg_bags = l_cabin / (PAX / 4) / 2
overhead_incr = mass_overhead / 12
x_cabin_start = x_front + len_front_storage
x_cg_over_front = (cg_start + x_cabin_start) / 2
longer_len_cg = cg_start + len_front_storage
def cg_shift(old_cg, old_mass, added_cg, added_mass):
    new_cg = (old_cg * old_mass + added_cg * added_mass) / (old_mass + added_mass)
    total_mass = old_mass + added_mass
    return new_cg, total_mass
front_fraction = (longer_len_cg - x_cabin_start) / 13
front_mass = front_fraction * mass_overhead
front_cargo_cg, front_cargo_m = cg_shift(x_cg_front, mass_front, x_cg_over_front, front_mass)
AC_cargo_f_cg, AC_cargo_f_m = cg_shift(longer_len_cg, oem, front_cargo_cg, front_cargo_m)
mass_aft = (1 - front_fraction) * mass_overhead
cg_aft = (x_cabin_start * 2 + 13) / 2
AC_cargo_aft_cg, AC_cargo_aft_m = cg_shift(longer_len_cg, oem, cg_aft, mass_aft)
AC_cargo_cg_1, AC_cargo_m_1 = cg_shift(AC_cargo_f_cg, AC_cargo_f_m, cg_aft, mass_aft)
AC_cargo_cg_2, AC_cargo_m_2 = cg_shift(AC_cargo_aft_cg, AC_cargo_aft_m, front_cargo_cg, front_cargo_m)
x_vals = [longer_len_cg, AC_cargo_f_cg, AC_cargo_aft_cg, AC_cargo_cg_1, AC_cargo_cg_2]
y_vals = [oem, AC_cargo_f_m, AC_cargo_aft_m, AC_cargo_m_1, AC_cargo_m_2]
loading_pax = 2 * m_pax / PAX
seat_incr = l_cabin / PAX * 4

def seating(where, start_cg, start_mass):
    cgs = start_cg#[x_vals[3]]
    masses = start_mass#[y_vals[3]]
    if where == "front":
        for i in range(0, 12):
            if i == 0:
                x_cg = x_cabin_start + 0.5 * seat_incr * (i + 0.5)
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            else:
                x_cg = x_cabin_start + seat_incr * (i + 0.5)
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            cgs.append(x_cg_new)
            masses.append(mass_new)
    elif where == "back":
        for i in range(0, 12):
            if i == 0:
                x_cg = x_cabin_start + (l_cabin - seat_incr * (i + 0.5))
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            else:
                x_cg = x_cabin_start + (l_cabin - seat_incr * (i + 0.5))
                x_cg_new, mass_new = cg_shift(cgs[i], masses[i], x_cg, loading_pax)
            cgs.append(x_cg_new)
            masses.append(mass_new)
    return cgs, masses

cgs_1, masses_1 = seating("front", [x_vals[3]], [y_vals[3]])
cgs_2, masses_2 = seating("back", [x_vals[3]], [y_vals[3]])
cgs_3, masses_3 = seating("front", [cgs_1[12]], [masses_1[12]])
cgs_4, masses_4 = seating("back", [cgs_1[12]], [masses_1[12]])

#LOAD FUEL
len_tank = 2.2
tank_mass = m_f #*2.4
cg_tank = x_cabin_start + l_cabin + 0.5 * len_tank
cg_end, mass_end = cg_shift(cgs_3[12], masses_3[12], cg_tank, tank_mass)

plt.figure()
plt.plot([x_vals[0], x_vals[1]], [y_vals[0], y_vals[1]], "black")
plt.plot([x_vals[0], x_vals[2]], [y_vals[0], y_vals[2]], "black")
plt.plot([x_vals[1], x_vals[3]], [y_vals[1], y_vals[3]], "black")
plt.plot([x_vals[2], x_vals[4]], [y_vals[2], y_vals[4]], "black")
plt.plot(cgs_1, masses_1, "red", marker = "s")
plt.plot(cgs_2, masses_2, "green", marker = "s")
plt.plot(cgs_3, masses_3, "blue", marker = "s")
plt.plot(cgs_4, masses_4, "orange", marker = "s")
plt.plot([cgs_3[12], cg_end], [masses_3[12], mass_end], "yellow", marker = "s")
plt.show()
print(cgs_1)
print(cgs_2)
print(masses_1)
print(masses_2)