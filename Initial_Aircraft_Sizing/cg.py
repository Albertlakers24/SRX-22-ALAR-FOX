import numpy as np
import scipy as sp
from Initial_Aircraft_Sizing.Fuselage import l_f, D_outer, S_f_wet, Mass_fuselage, l_nc, w_door_front, l_pax, w_lav, l_tank, prop_choice
from Initial_Aircraft_Sizing.Wing_planform import M_cruise, a_cruise, V_cruise, c_mac

# x_cg = Sum(mass * distance)/Sum(distance)
#Switches
#case = 1 # Case 1 = hydrogen combustion, 2 propellors on the wing, tank at aft
         # Case 2 = battery fuel hybrid in series/parrallel , 2 propellors on the fuselage, batteries at the center of the fuselage
         # Case 3 = battery fuel hybrid in series, 2 propellors at the wingtip, batteries at the center of the fuselage
         # Case 4 = hydrogen fuel cell, 2 propellors on the wing, tank at aft of fuselage

# Some parameters
M_Dive = M_cruise + 0.09     # Diving mach number
VDive = a_cruise * M_Dive
w_f = D_outer                # Width of fuselage
h_f = w_f                    # Height of fuselage, assumed to be a circular fuselage
x_cg_LEMAC = 0.25 * c_mac  # (predetermined value)
#distance of fuselage group components w.r.t. datum
x_fuselage = 0.4 * l_f
x_empennage = 0.9 * l_f
x_sys = 0.4 * l_f
x_wing = 0.4 * c_mac
x_prop_wing = -0.2 * c_mac
x_prop_fuse = 0.85 * l_f
x_hydrogen_tank = l_nc + w_door_front + l_pax + w_lav + l_tank/2
# l_batteries_center = 2.6 / np.pi*(2.5**2)                       # Size of batteries required assuming that the available cross sectional diameter for battery storage is 2.6 [m]
                                                                # both the cabin and the cargo compartment are filled with battery at the same lngitudinal position w.r.t.
                                                                # the aircraft
#l_batteries_center = 0.4 * l_f
l_batteries_wing = 0.75*c_mac-0.25*c_mac                        # assuming that if batteries are placed in wings, they will be placed in the space between the front (0.25c)
                                                                # and the rear spar (0.75c)
#x_batteries_center = l_nc + w_door_front + 0.5 * l_batteries_center     # cg position of the batteries if they're placed in the fuselage (fg)
x_batteries_center = 0.4 * l_f
x_batteries_wing = 0.5*c_mac                                            # cg position of the batteries if they're placed in the wings (wg)
x_batteries_forward = 0.3 * l_f
# Mass fraction of components w.r.t MTOM
if prop_choice == 1:           #Case 1 = hydrogen combustion, 2 propellors on the wing, tank at aft, more or less done
    Mf_fuselage = 0.14      # Mass fraction of components w.r.t MTOW
    Mf_prop = 0.09
    Mf_empennage = 0.04
    Mf_sys = 0.25
    Mf_wing = 0.13
    Mf_tank = 0.08   # Mass fraction of the tank to be updated

    print("======================================")
    print("For hydrogen combustion propulsion: ")
    print("hydrogen tank placement", x_hydrogen_tank/l_f)
    M_fg_frac = Mf_fuselage + Mf_empennage + Mf_sys + Mf_tank                                                   # sum of mass of the fuselage group
    M_wg_frac = Mf_wing + Mf_prop                                                                               # sum of mass of the wing group
    x_fg = (Mf_fuselage*x_fuselage+Mf_empennage*x_empennage+Mf_sys*x_sys+Mf_tank*x_hydrogen_tank) / M_fg_frac   # c.g. of the fuselage group
    x_wg_LEMAC = (Mf_wing * x_wing + Mf_prop * x_prop_wing) / M_wg_frac                                         # cg location of the wing group w.r.t. LEMAC
    x_LEMAC = x_fg - x_cg_LEMAC + M_wg_frac / M_fg_frac*(x_wg_LEMAC-x_cg_LEMAC)                                 # cg location of the lemac w.r.t. the  datum


elif prop_choice == 2:         # Case 2 = battery fuel hybrid in series/parrallel , 2 propellors on the fuselage, batteries at the center of the fuselage
    Mf_fuselage = 0.13   # Mass fraction of components w.r.t MTOM, need to be changed
    Mf_prop_fuse = 0.1
    Mf_prop_wing = 0.06
    Mf_sys = 0.25
    Mf_wing = 0.13
    Mf_empennage = 0.03
    Mf_batteries = 0.19
    print("======================================")
    print("For battery hybrid propulsion in series/ parallel configuration: ")
    M_fg_frac = Mf_fuselage + Mf_empennage + Mf_sys + Mf_batteries
    x_fg = ( Mf_fuselage * x_fuselage + Mf_empennage * x_empennage + Mf_sys * x_sys + Mf_batteries * x_batteries_center) / M_fg_frac
    M_wg_frac = Mf_wing + Mf_prop_wing
    x_wg_LEMAC = (Mf_wing * x_wing + Mf_prop_wing * x_prop_wing) / M_wg_frac                                                                 # cg location of the wing group w.r.t. to LEMAC
    x_LEMAC = x_fg - x_cg_LEMAC + M_wg_frac / M_fg_frac * (x_wg_LEMAC - x_cg_LEMAC)



elif prop_choice == 3:         # Case 3 = battery fuel hybrid in series, 2 propellors at the wingtip, batteries at the center of the fuselage
    Mf_fuselage = 0.12   ## Mass fraction of components w.r.t MTOM, need to be changed
    Mf_prop_fuse = 0.1
    Mf_prop_wing = 0.1
    Mf_sys = 0.24
    Mf_wing = 0.11
    Mf_empennage = 0.03
    Mf_batteries = 0.21
    print("======================================")
    print("For battery hybrid propulsion in series configuration: ")
    M_fg_frac = Mf_fuselage + Mf_empennage + Mf_sys + Mf_batteries
    x_fg = ( Mf_fuselage * x_fuselage + Mf_empennage * x_empennage + Mf_sys * x_sys + Mf_batteries * x_batteries_center) / M_fg_frac
    M_wg_frac = Mf_wing + Mf_prop_wing
    x_wg_LEMAC = (Mf_wing * x_wing + Mf_prop_wing * x_prop_wing) / M_wg_frac                                                                 # cg location of the wing group w.r.t. to LEMAC
    x_LEMAC = x_fg - x_cg_LEMAC + M_wg_frac / M_fg_frac*(x_wg_LEMAC-x_cg_LEMAC)                                 # cg location of the lemac w.r.t. the  datum



elif prop_choice == 4:          # Case 4 = hydrogen fuel cell, 2 propellors on the wing, tank at aft, , more or less done
    Mf_fuselage = 0.13  # Mass fraction of components w.r.t MTOM
    Mf_prop = 0.12
    Mf_empennage = 0.03
    Mf_sys = 0.25
    Mf_wing = 0.14
    Mf_tank = 0.05    # Mass fraction of the tank to be updated
    print("======================================")
    print("For hydrogen fuel cell architecture: ")
    print("hydrogen tank placement", x_hydrogen_tank/l_f)
    M_fg_frac = Mf_fuselage + Mf_empennage + Mf_sys + Mf_tank  # sum of mass of the fuselage group
    M_wg_frac = Mf_wing + Mf_prop  # sum of mass of the wing group
    x_fg = (Mf_fuselage * x_fuselage + Mf_empennage * x_empennage + Mf_sys * x_sys + Mf_tank * x_hydrogen_tank) / M_fg_frac  # c.g. of the fuselage group
    x_wg_LEMAC = (Mf_wing * x_wing + Mf_prop * x_prop_wing) / M_wg_frac
    x_LEMAC = x_fg - x_cg_LEMAC + M_wg_frac / M_fg_frac * (x_wg_LEMAC - x_cg_LEMAC) # The distance from nose to LEMAC.


x_oewcg = x_LEMAC + x_cg_LEMAC        # The distance from nose to cg.
print("The distance from zero point to LEMAC is ",x_LEMAC, "[m]")
print("The distance from zero point to cg is ", x_oewcg, "[m]")







