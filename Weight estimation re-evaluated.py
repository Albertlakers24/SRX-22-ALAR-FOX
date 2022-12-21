import numpy as np

# Outputs: MTOW OEW and WF
# Relationships: MTOW = WOE + WF + WPL
# Relationships: MTOW = WE + WF + WPLtot + Wtfo
# LOOK INTO MF/MTO FORMULA FROM ROLOEF READER
# LIQUID FUELS WILL HAVE THE SAME FORMULA AS JET-A FUEL
# BATTERIES WILL HAVE FORMULA WITHOUT ln()

# Manual inputs
Type = 1
range = 2

# Type 1: Hydrogen Combustion
# Type 2: Kerosene
# Range 1: 1,000 nmi (incl. contingencies)
# Range 2: 500 nmi (excl. contingencies)

# Handy unit conversion
lbs_kg = 0.453592  # 1 lbs = 0.453592 kg
ft_m = 0.3048  # 1 ft = 0.3048 m

# Constants
g = 9.80665  # gravity
f_con = 5 / 100  # -
Mres = 0.25  # -
e_kero = 46  # MJ/kg Specific Energy Kerosene
e_atj = 43.2  # MJ/kg Specific Energy SAF(ATJ)
e_lh2 = 120  # MJ/kg Specific Energy Liquid Hydrogen
e_bat = 1.1  # MJ/kg Specific Energy Battery (assuming 300Wh/kg)
eta_EM = 0.95  # - electric motor efficiency
eta_p_single = 0.80  # - propulsive efficiency - single engine/motor props
eta_p_twin = 0.82  # - propulsive efficiency - twin engine props
eta_p_turbo = 0.85  # - propulsive efficiency - regional turboprops
a_p = 0.5464  # - turboprop     WE = a*MTOW+b
a_j = 0.4985  # - turbojet
b_p = 1439 * g  # N turboprop
b_j = 1782.3 * g  # N turbojet
Psi = 0.0075  # - parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97  # - span efficiency factor (value based on Roelof reader p.46)
eta_eng_lh2 = 0.3  # Brake thermal efficiency 0.2-0.25
eta_eng_kero = 0.45  # Thermo efficiency
m_tank = 8.375  # kg
Number_tank = 68

# Constants depending on the aircraft
PAX = 50  # Number of passengers
V_cruise = 275 * 0.51444444  # m/s (TAS)  OUR:275      ATR42:300     ATR72: 280
h_cruise = 280 * 100 * ft_m  # m          OUR:280      ATR42: 7600m  ATR72: 200
R_norm = 1852000  # Range in m 1000 nmi
E = 30 * 60  # Loiter endurance in seconds
R_div = 185200  # m 100 nmi

# Constants depending on aircraft design
A = 12  # - our: 12 ATR72 value aspect ratio -> high A for slender wing
Cfe = 0.0045  # - (0.0045-0.005) equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1  # - (6.0-6.2) wetted area ratios -> depending on airframe structure

if Type == 1:
    e_f = e_lh2 * 1000000  # J/kg
    eta_p = eta_p_turbo  # -
    eta_eng = eta_eng_lh2  # -
    a = a_p  # - linear regression for OEW     (TYPES: turboprop, turbojet)
    b = b_p  # N linear regression for OEW     (TYPES: turboprop, turbojet)

    W_tanks = 1268 * 1.4 * g #m_tank * g * Number_tank  # N

if Type == 2:
    e_f = e_kero * 1000000  # J/kg
    eta_p = eta_p_turbo  # -
    eta_eng = eta_eng_kero  # -
    a = a_p  # - linear regression for OEW     (TYPES: turboprop, turbojet)
    b = b_p  # N linear regression for OEW     (TYPES: turboprop, turbojet)

# Calculations Aerodynamics
e = 1 / ((np.pi) * A * Psi + (1 / phi))  # -
Cd0 = Cfe * Swet_S  # -
CL = 0.72#np.sqrt(np.pi * Cd0 * A * e)  # -
CL_CD = 19.8
CD = 0.72 / CL_CD  # -

# Calculations Range
def Range(number):
    if number == 1:
        R_lost = 1 / 0.7 * (CL / CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm + R_lost) * (1 + f_con) + 1.2 * R_div + (E * V_cruise)  # m
    if number == 2:
        R_lost = 1 / 0.7 * (CL / CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm / 2 + R_lost) * (1 + f_con)  # m
    return R

# Calculations fuel fraction
mfuel_MTO_FULL = 1 - np.exp(-Range(range) / (eta_eng * eta_p * (e_f / g) * (CL / CD)))  # - normal
mbat_MTO_FULL = Range(range) / (eta_EM * eta_p * (e_bat / g) * (CL / CD))  # - battery
Mused = mfuel_MTO_FULL  # -

# Calculations Payload
WPAX = 200 * 0.453592 * PAX * g  # N
WPAXBAGGAGE = 40 * 0.453592 * PAX * g  # N Crew is bagageless
W_F_extra = 1.03  # 3%
WPLtot = WPAX + WPAXBAGGAGE  # N

if Type == 1:
    if range == 1:
        MTOW = (b + WPLtot + W_tanks) / (1 - a - (Mused * (1 + Mres)) * W_F_extra)
        WOE = a * MTOW + b + W_tanks
        WF = MTOW * (Mused * (1 + Mres)) * W_F_extra
        WF_extra = WF * (W_F_extra - 1)
    if range == 2:
        Mused = 1 - np.exp(-Range(1) / (eta_eng * eta_p * (e_f / g) * (CL / CD)))
        MTOW = (b + WPLtot + W_tanks) / (1 - a - (Mused * (1 + Mres)) * W_F_extra)
        WOE = a * MTOW + b + W_tanks
        Mused = 1 - np.exp(-Range(2) / (eta_eng * eta_p * (e_f / g) * (CL / CD)))
        WF = MTOW * (Mused)# * (1 + Mres))
        WF_extra = 0

if Type == 2:
    if Range == 1:
        MTOW = (b + WPLtot) / (1 - a - Mused * (1 + Mres))  # N
        WOE = a * MTOW + b  # N
        WF = MTOW * (Mused * (1 + Mres))  # N
    # if Range == 2:

#
# W_tanks =570*9.80665                                #N
#
# WF_extra = 38.04289576773986*9.81                   #N
#                #N
# print("Fuel at 500 nmi=", WF_halfrange/g, "kg")     #N

####ATR42
# WPLtot = WPAX+WPAXBAGGAGE                           #N
# WOE = 11750*g                                       #N
# MTOW = 18600*g                                      #N
# WF_halfrange = MTOW*(Mused*(1+Mres))                #N
# print("Fuel at 500 nmi=", WF_halfrange/g, "kg")     #N


print("MTOW =", MTOW, "N")
print("Fuel weight =", WF, "N")
print("Operational Empty Weight =", WOE, "N")
print("WPLtot =", WPLtot, "N")
print("WF_extra =", WF_extra, "N")
print("---------------------------------------")
print("m_TO=", MTOW / g, "kg")
print("Fuel mass=", WF / g, "kg")
print("m_OE=", WOE / g, "kg")
print("m_PLtot=", WPLtot / g, "kg")
print("m_tanks =", W_tanks / g, "kg")
print("mF_extra =", WF_extra / g, "kg")

print("efficiency =", eta_eng * eta_p)

print("CL/CD =", CL / CD)
print("Cl=", CL)
print("CD", CD)
