import numpy as np
import scipy as sp
from Class_I_Weight_Estimation import b, A, Swet_S, MTOW
from Lift_Drag_Polar import p

gamma = 1.4
taper = 0.4                            # Taper ratio (from 0 to 1), 0.4 for unswept wing. check Raymer
c_r = 2*Swet_S/((1+taper)*b)           # Tip chord
c_t = c_r * taper                      # Root chord
qhat = 0.5 * gamma * p * Mcruise             # Cruise Mach Speed
C_Lhat = MTOW/(qhat*Swet_S)            # Cruise lift coefficient
t_c_ratio = min(0.18, C_Lhat)
