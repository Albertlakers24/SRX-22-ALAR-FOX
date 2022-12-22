import numpy as np
from Class_I_Weight_Estimation.Drag_Polar_new import max_CL_CD
from Constants import *
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *

# Outputs: MTOW OEW and WF
# Relationships: MTOW = WOE + WF + WPL
# Relationships: MTOW = WE + WF + WPLtot + Wtfo
# LOOK INTO MF/MTO FORMULA FROM ROLOEF READER
# LIQUID FUELS WILL HAVE THE SAME FORMULA AS JET-A FUEL


#Constants
f_con = 5/100                   #-
m_res = 0.25                    #-
a_regression = 0.5088           #- Regression value
b_regression = 1199.7           #- Regression value
CL_CD = round(max_CL_CD)        #- Imported value
total_eff = eta_fuelcell * eta_inverter * eta_wire * eta_EM * eta_prop #one inverter
m_f_extra = 0.03                #-
fuel_mass_ref = 669             #kg
m_tanks = 1.4 * fuel_mass_ref   #kg
fc_power_density = 3            #kW/kg
inverter_power_density = 30     #kW/kg
em_power_density = 15           #kW/kg

#Calculations Range
def mf_mMTO(range):
    if range == "full":
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm + R_lost) * (1 + f_con) + 1.2 * R_div + (E * V_cruise)  # m
    elif range == 500:
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm / 2 + R_lost) * (1 + f_con)  # m
    mf_mMTO_fraction = 1 - np.exp((-1 * round(R)) / (total_eff * (e_lh2 / g) * CL_CD))
    return mf_mMTO_fraction

# print(mf_mMTO("full"))
m_payload = m_pax + m_baggage
m_mto = (b_regression + m_payload) / (1 - a_regression - ((mf_mMTO("full") * 2.4) * (1 + m_res)) * (1 + m_f_extra))
m_oem = a_regression * m_mto + b_regression

# e = 2
# while e > 1.037:
#     power_shaft = W_P_design * m_mto
#     fc_mass = (power_shaft / total_eff * eta_fuelcell * eta_prop) / fc_power_density
#     em_mass = (power_shaft / eta_EM) / em_power_density
#     inverter_mass = (power_shaft / eta_EM / eta_wire / eta_inverter) / inverter_power_density
#     m_mto_old = m_mto
#     m_mto = (fc_mass + em_mass + inverter_mass) + m_mto_old
#     print(m_mto, m_mto_old)
#     e = (m_mto / m_mto_old)

# print(m_mto, "m_mto")



def fuel_mass(range):
    if range == "full":
        m_f = m_mto * (mf_mMTO(range) * (1 + m_res)) * (1 + m_f_extra)
    elif range == 500:
        m_f = m_mto * (mf_mMTO(range)) * (1 + m_f_extra)
    return m_f

# print(fuel_mass("full"))
# print(fuel_mass_ref)
#
power_needed = m_mto * g / W_P_design / 10**3
print(power_needed, m_mto)




