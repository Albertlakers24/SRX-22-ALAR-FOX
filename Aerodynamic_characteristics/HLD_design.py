import numpy as np
from Constants import *
from Wing_lift_estimation import Calculate_wingsweep
from Initial_Aircraft_Sizing.Wing_planform import Sw,b,c_mac,c_r,taper
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
Cl_max_airfoil = 1.61
CL_over_Cl_ratio = 0.9

#Design Options
te_hld = ["single slotted","double slotted","fowler"]
le_hld = ["None","LE_flap","slat"]
lm = 1.85                                               #Minimum distance from center line to TE HLD [m]
lm_LE = 4                                               #Minimum distance from center line to LE HLD [m]

#Planform data
LE_hinge_line_angle_deg = Calculate_wingsweep(0,0.15,taper)   #front spar position at 15% chord
TE_hinge_line_angle_deg = Calculate_wingsweep(0,0.6,taper)    #rear spar position at 60% chord
LE_angle_deg = Calculate_wingsweep(0,0,taper)                 #Leading edge angle
TE_angle_deg = Calculate_wingsweep(0,1,taper)                 #Trailing edge angle
aileron_percent = 0.6                                   #aileron starting position

file = open("HLD_design_choice.txt","w")
alpha_angle = np.radians(90 - TE_angle_deg)
beta_angle = np.radians(90+LE_angle_deg)

A1 = trapezoid_area(lm,l_top(alpha_angle,beta_angle,lm,c_r),c_r)
lm_top = l_top(alpha_angle,beta_angle,lm,c_r)
l3_sequence = yield_func(lm_LE,(b/2)*aileron_percent,0.05)
try:
    file.write("ld[m]\tld[%]\tSwfS\tCf\tDCL\tDCLalpha\tdf\tc_prime_over_c\tflap_type\n")
    deltaCL_max_corrected = CL_max_req - Cl_max_airfoil * CL_over_Cl_ratio
    A1 = trapezoid_area(lm,l_top(alpha_angle,beta_angle,lm,c_r),c_r)
    lm_top = l_top(alpha_angle,beta_angle,lm,c_r)
    l3_sequence = yield_func(lm_LE,(b/2)*aileron_percent,0.05)
    for l3 in l3_sequence:
        A2 = trapezoid_area(l3,l_top(alpha_angle,beta_angle,l3,lm_top),lm_top)
        if 0.8 <= (A2-A1)/(Sw/2):
            ld = l3
            ld_percent = l3/(b/2) * 100
            break
        else:
            ld = -1
            ld_percent = -1
    print(lm_LE, " ", ld)
    for flap_type in te_hld:
        for LE_type in le_hld:
            deltaf = yield_func(30,45,5)
            print(flap_type)
            for df in deltaf:
                SwfS_sequence = yield_func(0.3,0.9,0.01)
                for SwfS in SwfS_sequence:
                    Cf_sequnce = yield_func(0.10,0.41,0.01)
                    for Cf in Cf_sequnce:
                        DPS = False
                        c_prime_over_c_LE_sequence = yield_func(1, 1.1, 0.01)
                        for c_prime_over_c_LE in c_prime_over_c_LE_sequence:
                            '''Start of looped code'''
                            deltaClmax_HLD, c_prime_over_c = HLD_TE_deltaClmax(Cf, df,
                                                                                  flap_type)  # Calculating delta Cl max for trailing edge HLD
                            deltaClmax_LE = HLD_LE_deltaClmax(LE_type,
                                                                c_prime_over_c_LE)  # Calculating delta Cl max for leading edge HLD

                            DP = 0.9 * deltaClmax_HLD * SwfS * np.cos(np.radians(TE_hinge_line_angle_deg)) + 0.9 * deltaClmax_LE * 0.75 * np.cos(
                                np.radians(LE_hinge_line_angle_deg))  # Calculating total CL max increase
                            DPA = -15 * SwfS * np.cos(np.radians(LE_hinge_line_angle_deg))  # Calculating change in CL alpha
                            if DP >= deltaCL_max_corrected:
                                l1_sequence = yield_func(lm,b*aileron_percent,0.05)
                                for l1 in l1_sequence:
                                    A2 = trapezoid_area(l1, l_top(alpha_angle, beta_angle, l1, lm_top), lm_top)
                                    if SwfS <= (A2-A1)/ (Sw/2):
                                        ld = l1
                                        ld_percent = l1/(b/2)*100
                                        break
                                    else:
                                        ld = -1
                                        ld_percent = -1
                                file.write(str(round(ld,2)) + '\t' +str(round(ld_percent,2)) + '\t' +str(round(SwfS,2)) + '\t' + str(round(Cf,2))+ '\t' + str(round(DP,2))+ '\t' + str(round(DPA,2))+ '\t' + str(df)+ '\t' + str(round(c_prime_over_c,2))+ '\t' + flap_type + '+' + LE_type + '\n')
                                DPS = True
                            if DPS:    break
                        if DPS:     break
finally:
    file.close()
    print("end code")