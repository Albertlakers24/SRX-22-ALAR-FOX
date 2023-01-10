import numpy as np
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import V_approach_stall
from Initial_Aircraft_Sizing.Wing_planform import Sw
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

#def HLD_TE_deltaClmax(Cf, df, flap_type):
CL_max_land_design = 1.1 * CL_max_landing
print(CL_max_land_design)