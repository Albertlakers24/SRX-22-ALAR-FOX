import numpy as np
from Constants import *
from Wing_lift_estimation import Calculate_wingsweep
from Initial_Aircraft_Sizing.Wing_planform import Sw,b,c_mac,c_r
def yield_func(low,high,step):
    i = low
    while i <= high:
        yield i
        i+=step
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

def trapezoid_area(h,top,base):
    area = (top+base)/2 * h
    return area

def l_top(alpha,beta,height,base):
    x1 = height/ np.tan(alpha)
    x2 = height/ np.tan(beta)
    top = base - x1 - x2
    return top

#Lift Data
CL_max_req = 1.1 * CL_max_landing
Cl_max_airfoil = 1.5
CL_over_Cl_ratio = 0.88

#Design Options
flaps = ["single slotted","double slotted","fowler"]
leading = ["None","LE_flap","slat"]
lm = 1.85                                               #Minimum distance from center line to TE HLD [m]
lm_LE = 4                                               #Minimum distance from center line to LE HLD [m]

#Planform data
LE_hinge_line_angle_deg = Calculate_wingsweep(0,0.15)   #front spar position at 15% chord
TE_hinge_line_angle_deg = Calculate_wingsweep(0,0.6)    #rear spar position at 60% chord
LE_angle_deg = Calculate_wingsweep(0,0)                 #Leading edge angle
TE_angle_deg = Calculate_wingsweep(0,1)                 #Trailing edge angle
aileron_percent = 0.6                                   #aileron starting position

#file = open("HLD_design_choice.txt","w")
#file.write("ld[m]\tld[%]\tSwfS\tCf\tDCL\tDCLalpha\tdf\tc_prime_over_c\tflap_type\n")
#file.close()
deltaCL_max_corrected = CL_max_req - Cl_max_airfoil*CL_over_Cl_ratio
alpha_angle = np.radians(90 - TE_angle_deg)
beta_angle = np.radians(90+LE_angle_deg)

A1 = trapezoid_area(lm,l_top(alpha_angle,beta_angle,lm,c_r),c_r)
lm_top = l_top(alpha_angle,beta_angle,lm,c_r)
l3_sequence = yield_func(lm_LE,(b/2)*aileron_percent,0.05)
for l3 in l3_sequence:
    A2 = trapezoid_area(l3,l_top(alpha_angle,beta_angle,l3,lm_top),lm_top)
    if 0.75 <= (A2-A1)/(Sw/2):
        ld = l3
        ld_percent = l3/b * 100
        break
    else:
        ld = -1
        ld_percent = -1
print(lm_LE, " ", ld)