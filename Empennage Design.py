import numpy as np
from Class_I_Weight_Estimation import V_cruise, MTOW
from Lift_Drag_Polar import b
from Wing_planform import Sw, c_mac, specific_gas_constant, T, p, gamma, M_cross

#Constants from other files
cw = c_mac                          # m c_mac from Wing_planform
bw = b                              # m b from Lift_Drag_Polar
gamma = 1.4                         # Specific heat ratio of gas
M_cross = 0.935                     # Technology factor for supercritical airfoils

#CHOOSE
lh = 24                             # m from CG calculations I think
lv = 24                             # m from CG calculations
Vh = 0.68113                        # TBD horizontal tail volume from data  -> 0.68113 (propeller)
Vv = 0.08                           # TBD vertical tail volume from data    -> 0.08 (expected for propeller)
sweeph_quarter = 0                  # degrees (0 for propeller)
sweepv_LE = 25                      # degrees (0-50) propeller
Av = 1.3                            # - (1-2)
Ah = 4                              # - (3-5)
taperv = 0.5                        # - (0.3-0.7) taper ratio
taperh = 0.45                       # - (0.3-1) taper ratio
sweepv_half = 25                    # degrees -> depending on the t/c
fraction_df = 0.19                  # -
sweep_df_LE = 74                    # degrees TBD

#Check with graph
Sh_stability = 1.06                 # Determine from the stability graph

#Intermediate Equations --now for the wing, to be changed?!!!!!!!
l_opt =


####HORIZONTAL TAIL
#Calculations changes for Conventional =1, Conventional plus dorsal =2
Type = 1

if Type ==1:
    Sh = Vh * (Sw * cw) / lh                # m^2 Horizontal tail surface area
    if abs((Sh-Sh_stability)/Sh)*100 < 10:
        #Optimum Tail moment arm

        #Determine span, root chord, tipchord MAC
        bh = np.sqrt(Sh*Ah)                 # m horizontal tail span
        c_rh = (2*Sh)/((1+taperh)*bh)       # m root chord
        c_th = c_rh * taperh                # m Root chord
        t_c_ratioh = 0.18                   #ASSUMPTION min(0.18, ((M_cross - M_dd) - 0.115 * (C_Lhat ** 1.5)))        # no sweep in horizontal tail! thickness to chord ratio
        c_mach_h = (2 / 3) * c_rh * ((1 + taperh + taperh ** 2) / (1 + taperh))     # length of MAC
        y_mach_h = 0.5 * (1 / 3) * (1 + 2 * taperh) / (1 + taperh) * bh             # Spanwise location of MAC
        #Position the surface

        #Print Statements
        print("Sh =", Sh, "m^2")
        print("")
    else:
        print("Sh =", Sh, "m^2")
        print("Rerun Code for lh")

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
    print("H")''''