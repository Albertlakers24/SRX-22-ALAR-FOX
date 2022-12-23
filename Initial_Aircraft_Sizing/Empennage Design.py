import numpy as np
from Wing_planform import Sw, c_mac, specific_gas_constant, gamma, M_cross,b



##Design 1; Fuel Cell
##Design 2; Liquid Hydrogen Combustion
##Design 3; Hybrid Electric Parallel Series
##Design 4; Hybrid Electric Series CHeck
Design = 3

#Constants from other files
D_outer = 2.8 + 0.21*2              #m
D_inner = 2.8                       #m

if Design ==1:                  #FUEL CELL -> DONE
    ##AIRCRAFT DESIGN fixed inputs
    cw = 2.2                # from Gabriel (double tapered)
    bw = 25                 # from Gabriel
    Sw = 53                 # from Gabriel
    l_f = 23.876            # m
    xcg_aft = 11.28         # m
    l_ultimate = 19.755      # m

    ##Change
    x_h = 22.3              #m
    x_v = 20.9              #m
    Vh = 0.86               # TBD horizontal tail volume from data  -> 0.68113 (propeller)
    Vv = 0.08               # TBD vertical tail volume from data    -> 0.08 (expected for propeller)
    Av = 1.3                # - (1-2)
    Ah = 4                  # - (3-5)
    taperv = 0.5            # - (0.3-0.7) taper ratio
    taperh = 0.75           # - (0.3-1) taper ratio
    Kc = 1.4

    l_v = x_v - xcg_aft     #
    l_h = x_h - xcg_aft     # m
    l_opt = Kc * np.sqrt((4 * c_mac * Sw * Vh) / (np.pi * D_outer))

if Design ==2:              #Liquid hydrogen combustion
    #AIRCRAFT DESIGN FROM GABRIEL
    cw = 2.3                    # m
    bw = 27                     # m
    Sw = 60                     # m^2
    l_f = 26.8                  # m
    xcg_aft = 12.79             # m
    l_ultimate = 22.679         # m   - update

    Vh = 0.9                    # TBD horizontal tail volume from data  -> 0.68113 (propeller)
    Vv = 0.068                  # TBD vertical tail volume from data    -> 0.08 (expected for propeller)
    Av = 1.3                    # - (1-2)   little lower because less propellers on the wing
    Ah = 4.5                    # - (3-5)
    taperv = 0.5                # - (0.3-0.7) taper ratio
    taperh = 0.6                # - (0.3-1) taper ratio
    x_v = 23.9
    x_h = x_v +1.5

    l_h = x_h-xcg_aft
    l_v = x_v-xcg_aft

    Kc = 1.2
    print("-----------Location check---------")
    print("l_v =", l_v)
    l_opt = Kc * np.sqrt((4 * c_mac * Sw * Vh) / (np.pi * D_outer))
    print("l_opt=", l_opt)

if Design ==3:                          #Parallel series hybrid
    # AIRCRAFT DESIGN FROM GABRIEL
    cw = 2.9                    # m
    bw =  33                    # m
    Sw =   92                   # m^2
    l_f = 21.7                  # m
    xcg_aft = 9.44              # m
    l_ultimate = 16.96          # m-----update

    Vh = 0.9                    # TBD horizontal tail volume from data  -> 0.68113 (propeller)
    Vv = 0.058                   # TBD vertical tail volume from data    -> 0.08 (expected for propeller)
    Av = 1.5                    # - (1-2) ATR42 ->1.5
    Ah = 4.5                    # - (3-5) ATR42 ->4.5
    taperv = 0.4                # - (0.3-0.7) taper ratio
    taperh = 0.6                # - (0.3-1) taper ratio
    x_v = 17.55                  # m
    x_h = x_v+2.4               # m

    l_h = x_h - xcg_aft
    l_v = x_v - xcg_aft

    Kc = 1.2
    print("-----------Location check---------")
    print("l_v =", l_v)
    l_opt = Kc * np.sqrt((4 * c_mac * Sw * Vh) / (np.pi * D_outer))
    print("l_opt=", l_opt)

if Design ==4:                  #Hybrid; series
    # AIRCRAFT DESIGN FROM GABRIEL
    cw =  3.1              # m
    bw =  36                # m
    Sw = 106               # m^2
    l_f =21.7             # m
    xcg_aft = 9.32         # m
    l_ultimate = 14         # m  -----update

    Vh = 0.9                # TBD horizontal tail volume from data  -> 0.68113 (propeller)
    Vv = 0.06               # TBD vertical tail volume from data    -> 0.08 (expected for propeller)
    Av = 1.5                # - (1-2)
    Ah = 4.5                # - (3-5)
    taperv = 0.5            # - (0.3-0.7) taper ratio
    taperh = 0.8            # - (0.3-1) taper ratio
    x_v = 18                # m
    x_h = x_v+2.5           # m

    l_h = x_h - xcg_aft
    l_v = x_v - xcg_aft

    Kc = 1.2
    print("-----------Location check---------")
    print("l_v =", l_v)
    l_opt = Kc * np.sqrt((4 * c_mac * Sw * Vh) / (np.pi * D_outer))
    print("l_opt=", l_opt)

#Intermediate Equations
print("l_f=", l_f)

####HORIZONTAL TAIL
#Calculations changes for Conventional =1, Conventional plus dorsal =2
switch =1
Type =1

if Type ==1:
    Sh = Vh * (Sw * cw) / l_h                # m^2 Horizontal tail surface area
    if switch ==1:         #abs((Sh-Sh_stability)/Sh)*100 < 10:
        #Determine span, root chord, tipchord MAC
        bh = np.sqrt(Sh*Ah)                 # m horizontal tail span
        c_rh = (2*Sh)/((1+taperh)*bh)       # m root chord
        c_th = c_rh * taperh                # m Root chord
        t_c_ratioh = 0.18                   #ASSUMPTION min(0.18, ((M_cross - M_dd) - 0.115 * (C_Lhat ** 1.5)))        # no sweep in horizontal tail! thickness to chord ratio
        c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))     # length of MAC
        y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh             # Spanwise location of MAC
        #Position the surface

        #Print Statements
        print("-----------Horizontal Tail Sizing------------")
        print("Sh =", Sh, "m^2")
        print("bh =", bh, "m")
        print("cr_h =", c_rh, "m")
        print("ct_h =",c_th, "m")
        print("c_mach =", c_mach_h)
        print("y_machh =", y_mach_h)
        print("ratio wing areas =", Sh/Sw)
        print("tail_arm_h=", x_h - xcg_aft)

        print("Av=", Av)
        Sv = Vv * (Sw * bw) / l_v                           #correct
        bv = np.sqrt(Av*Sv)                                 #correct
        c_mac_v = Sv/bv                                     #correct
        c_rv = 2*Sv/((1+taperv)*bv)
        c_tv = taperv*c_rv

        #c_rv = (2+Sv)/((1+taperv)*bv)       #

        #c_mac_v =(2/3)*c_rv*((1+taperv+taperv**2)/(taperv+1))
        y_mach_v = 0.5 * (1 / 3) * (1 + 2 * taperv) / (1 + taperv) * bv
        print("-----------Vertical Tail Sizing------------")
        print("Sv =", Sv, "m^2")
        print("bv =", bv, "m")
        print("cr_v =", c_rv, "m")
        print("ct_v =",c_tv, "m")
        print("c_macv =", c_mac_v)
        print("y_machv =", y_mach_v)
        print("tail_arm_v=", x_v-xcg_aft)

        print("Sh ratio", Sw/Sh)
        print("Sv ratio", Sw/Sv)

        print("----------------DRAWING CALCULATIONS-------------")
        #print("Location Horizontal Tail =", xcg_aft + l_opth)
        #print("Distance from TE to tail =", l_f - (xcg_aft + l_opth))
        print("l_h", l_h)
        print("l_v", l_v)
        print("part outside", (x_v+(3/4)*c_rv)-l_f)
        print("x_h=", x_h)
        print("x_v =", x_v)
        print("nose to LE vertical tail=", x_v-(1/4*c_rv))
        print("nose to LE horizontal tail=", x_h-(1/4*c_rh))
        if (x_v-(1/4)*c_rv)-l_ultimate < 0:
            print("we are interfering :( by", (x_v - (1 / 4) * c_rv) - l_ultimate)
        else:
            print("not interfering by:",(x_v - (1 / 4) * c_rv) - l_ultimate)
    else:
        print("Sh =", Sh, "m^2")
        print("Rerun Code for lh")



#print("distance outside fuselage=", l_f-(xcg+l_opt+c_rh))
#Check with graph
#Sh_stability = 11.7                 # Determine from the stability graph
#CHOOSE
#lh = 24                            # m from CG calculations I think
#lv = 24                            # m from CG calculations
#sweep_df_LE = 74                    # degrees TBD
#sweeph_quarter = 0                  # degrees (0 for propeller)
#sweepv_LE = 25                      # degrees (0-50) propeller
#sweepv_half = 25                    # degrees -> depending on the t/c
#fraction_df = 0.19                  # -






''''
if Type ==2:
    Sh = Vh * (Sw * cw) / lh  # m^2 Horizontal tail surface area
    if abs((Sh - Sh_stability) / Sh) * 100 < 10:
        # Determine span, root chord, tipchord MAC
        bh = np.sqrt(Sh * Ah)  # m horizontal tail span
        c_rh = (2 * Sh) / ((1 + taperh) * bh)  # m root chord
        c_th = c_rh * taperh  # m Root chord
        t_c_ratioh = min(0.18, ((M_cross - M_dd) - 0.115 * (
                    C_Lhat ** 1.5)))  # no sweep in horizontal tail! thickness to chord ratio
        c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))  # length of MAC
        y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh  # Spanwise location of MAC
        #Size Dorsal
        S_df  ##Calculate!!!
        h_df = np.sqrt((2 * S_df) / (np.tan(sweep_df_LE) - np.tan(sweepv_LE)))


        # Position the surface

        # Print Statements
        print("Sh =", Sh, "m^2")

if Type ==3:
    print("H")


#VERTICAL TAIL

if Type ==1:
    #Choose longitudinal location of vertical tail aerodynamic center
    #determine moment arm to the most aft CG position - lv
    #Compute required tail area:
    Sv = Vv*(Sw*bw)/lv                      # m^2 Horizontal tail surface area
    #Determine sweep, AR, taper
    #Determine span, cr, ct, MAC
    bv = np.sqrt(Sv * Av)  # m horizontal tail span
    c_rv = (2 * Sv) / ((1 + taperv) * bv)  # m root chord
    c_tv = c_rv * taperv  # m Root chord
    t_c_ratiov = min(0.18, (((np.cos(sweepv_half))**3)*(M_cross - M_dd) - 0.115 * (C_Lhat ** 1.5))/(np.cos(sweepv_half))**2)  # no sweep in vertical tail! thickness to chord ratio
    c_machv = (2 / 3) * c_rv * ((1 + taperv + taperv ** 2) / (1 + taperv))  # length of MAC
    y_machv = 0.5 * (1 / 3) * (1 + 2 * taperv) / (1 + taperv) * bv
    #position the surface

if Type ==2:
    print("goodbye")
if Type ==3:
    print("H")'''