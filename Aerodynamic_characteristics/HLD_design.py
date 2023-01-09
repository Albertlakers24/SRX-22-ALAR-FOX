import numpy as np
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Initial_Aircraft_Sizing.Wing_planform import Sw
def HLD_TE_deltaClmax(Cf,df,flap_type):
    if flap_type == "single slotted":
        delta_c =       #Torenbeek page 533
        c_prime_over_c = 1 + delta_c*Cf
        delta_clmax = 1.3
    if flap_type == "double slotted":

    if flap_type == "fowler":

def HLD_TE_deltaClmax(Cf, df, flap_type):