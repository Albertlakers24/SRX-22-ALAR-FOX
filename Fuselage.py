import numpy as np
import scipy as sp

# Some constants
k_cabin = 1.08      #For single aisle

# Dimensions under the regulation
w_door_front = 0.61
w_door_back = 0.508
l_lav = 0.914
w_lav = 0.914
l_galley = 0.762
w_galley = 0.914

# Design choices
n_SA = 3            # Number of seats abreast
n_PAX = 48          # Number of passengers (46-50)
n_row = np.ceil(n_PAX/n_SA)





# Outputs
D_eff = np.sqrt(width*height)

