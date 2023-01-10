import numpy as np
from Constants import *
from Wing_lift_estimation import Calculate_wingsweep
def HLD_TE_deltaClmax(Cf,df,flap_type):
    if flap_type == "single slotted":
        delta_c = 0      #Torenbeek page 533
        c_prime_over_c = 1 + delta_c*Cf
        delta_clmax_TE = 1.3
    if flap_type == "double slotted":
        if df < 15.01:
            delta_c = 0.3/15                                                    #Torenbeek page 533
        else:
            delta_c = 0.3 + df * ((0.70 - 0.30) / (60 - 15))                    #Torenbeek page 533
        c_prime_over_c = 1 + delta_c * Cf
        delta_clmax_TE = 1.3 * c_prime_over_c
    if flap_type == "fowler":
        if df < 10.01:
            delta_c = 0.45 / (1/1.2 * 10)                                       #Torenbeek page 533
        else:
            delta_c = 0.45 + df * ((0.65 - 0.45 )/(45 - (1/1.2*10)))            #Torenbeek page 533
        c_prime_over_c = 1 + delta_c*Cf
        delta_clmax_TE = 1.3 * c_prime_over_c

    return delta_clmax_TE,c_prime_over_c

def HLD_LE_deltaClmax(flap_type,c_prime_over_c):
    if flap_type == "LE_flap":
        delta_clmax_LE = 0.3
    if flap_type == "slat":
        delta_clmax_LE = 0.4 * c_prime_over_c
    else:
        delta_clmax_LE = 0
    return delta_clmax_LE

#Lift Data
CL_max_req = 1.1 * CL_max_landing
Cl_max_airfoil = 1.2
CL_over_Cl_ratio = 0.88

#Design Options
flaps = ["single slotted","double slotted","fowler"]
leading = ["None","LE_flap","slat"]

