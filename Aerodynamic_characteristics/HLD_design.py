import numpy as np
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Initial_Aircraft_Sizing.Wing_planform import Sw
CL_max_required = (m_mto * g) / (1/2 * ISA_calculator(takeoff_critical,0)[2]*Sw*V_stall)