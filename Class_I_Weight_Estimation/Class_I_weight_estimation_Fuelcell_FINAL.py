import numpy as np
from Class_I_Weight_Estimation.Drag_Polar_new import max_CL_CD
from Class_I_Weight_Estimation.Wing_Loading_Diagram import *

#Constants
f_con = 5/100                   #-
m_res = 0.25                    #-
a_regression = 0.5617           #- Regression value
b_regression = 1199.7           #- Regression value
CL_CD = round(max_CL_CD)        #- Imported value
total_eff = eta_fuelcell * eta_inverter * eta_wire * eta_EM * eta_prop #one inverter
m_f_extra = 0.03                #-
fuel_mass_ref = 669             #kg
m_tanks = 1.4 * fuel_mass_ref   #kg

#Calculations Range
def mf_mMTO(range):
    if range == "full":
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm + R_lost) * (1 + f_con) + 1.2 * R_div + (E * V_cruise)  # m
    elif range == 500:
        R_lost = 1 / 0.7 * (CL_CD) * (h_cruise + (V_cruise ** 2 / (2 * g)))  # m
        R = (R_norm / 2 + R_lost) * (1 + f_con)  # m
    mf_mMTO_fraction = 1 - np.exp((-1 * R) / (total_eff * (e_lh2 / g) * CL_CD))
    return mf_mMTO_fraction

m_payload = m_pax + m_pax_baggage
m_crew_total = m_crew + m_crew_baggage

fc_mass = (1 / total_eff * eta_fuelcell * eta_prop) / fc_power_density
em_mass = (1 / eta_EM) / em_power_density
inverter_mass = (1 / eta_EM / eta_wire / eta_inverter) / inverter_power_density
masses_sum = (fc_mass + em_mass + inverter_mass) / W_P_design * g / 10**3

def mtom(oew_ratio, range):
    mtom = (m_payload + m_crew_total) / (1 - 2.4 * mf_mMTO(range) * (1 + m_res) - oew_ratio)
    return mtom
oew_mtom = a_regression + masses_sum
m_mto = mtom(oew_mtom, "full")

def fuel_mass(oew_ratio, range):
    if range == "full":
        m_f = mtom(oew_ratio, range) * (mf_mMTO(range) * (1 + m_res)) * (1 + m_f_extra)
    elif range == 500:
        m_f = mtom(oew_ratio, range) * (mf_mMTO(range)) * (1 + m_f_extra)
    return m_f
m_f = fuel_mass(oew_mtom, "full")
oem = oew_mtom * m_mto + m_f * 1.4
print("MTOM:", m_mto)
print("OEM:", oem)
print("Fuel mass:", m_f)