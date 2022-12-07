import numpy as np
import scipy as sp
#Switch
switch = 1          # put 1 for circular fuselage, 2 for double bubble fuselage

# Some constants
k_cabin = 1.08      #For single aisle

# Dimensions under the regulation
w_door_front = 0.61
w_door_back = 0.508
l_lav = 0.914
w_lav = 0.914
l_galley = 0.762
w_galley = 0.914
seat_pitch = 1

# Design choices
n_SA = 4           # Number of seats abreast

if n_SA == 3:
    n_PAX = 48                  # Number of passengers (46-50)
    D_inner = 2.48              # inner diameter of the fuselage cross section [m]
elif n_SA == 4:
    n_PAX = 48                  # Number of passengers (46-50)
    D_inner = 2.8               # inner diameter of the fuselage cross section [m]
elif n_SA == 5:
    n_PAX = 50                  # Number of passengers (46-50)
    D_inner = 3.5               # inner diameter of the fuselage cross section [m]
n_row = np.ceil(n_PAX/n_SA)     # number of rows
l_cp = 4                        # cockpit length [m]
double_bubble_height = 2.7

if switch == 1:
    height = D_inner
else:
    height = double_bubble_height

# Basic Length/Dimensions calculations
D_eff = np.sqrt(height*D_inner)             # effective (inner) diameter [m]
l_t = 1.6 * D_eff                           # Length of tail [m]
l_cabin = n_row * k_cabin                   # Length of the cabin [m]
skin_thickness = 0.084 + 0.045 * D_eff      # fuselage skin thickness [m]
D_outer = D_eff + skin_thickness            # outer fuselage diameter [m]

fuel_at_belly = True

if fuel_at_belly == True:
    l_tank = 0
else:
    l_tank = 4                                      # Length of the fuel tank as fuselage section
l_f = l_cabin + l_t + l_cp + l_tank                 # length of the fuselage [m]
fineness_ratio = l_f/D_eff                          # Fineness ratio
l_nc = 2 * D_eff                                    # length of the nose cone [m]
l_tc = 3 * D_eff                                    # length of the tail cone
# Skin Drag Computation



# Fuselage mass computation



#printing results
print("N_row:", n_row)
print("N_SA:", n_SA)
print("D_in [m]:", D_inner)
print("D_out [m]:", D_outer)
print("skin thickness [m]:", skin_thickness)
print("L_tank [m]:", l_tank)
print("l_fuselage [m]:", l_f)
print("Fineness Ratio: ", fineness_ratio)

