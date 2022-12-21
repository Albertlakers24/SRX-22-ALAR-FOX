import numpy as np
from matplotlib import pyplot as plt


g = 9.80665

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 85200   # in m ---> 100nmi
t_E =  45 * 60    # in seconds - endurance time
m_pldes = 5443  # Payload mass designed for

# Propulsion charecteristics
eta_p = 0.85                             #Propulsive efficiency twin turboprop
e_f = 120 * 10**6

# Masses Aircraft --> Max, empty, fuel - per propulsion system
def configuration_values(prop_type):
    if prop_type == 1:  # LH2 fuel cell
        m_f = 669 # LH2 mass in kg
        m_tank = 936 # 1.4*m_fh
        m_oe = 12195 + m_tank# Operating empty mass + tank mass
        m_mto = 19243 #19370
        pl_increase = 1.06145
        LD_crs = 22  # Lift to drag ratio during cruise
        eta_eng = 0.6 * 0.97 * 0.995 * 0.85 * 0.95
    if prop_type == 2:    # LH2 Combustion
        m_mto = 21728 # in kg --> Max take off
        m_oe = 15086  # in kg --> Operating empty
        m_f = 1268  # in kg --> Fuel mass
        pl_increase = 1.116
        LD_crs = 19.8  # Lift to drag ratio during cruise
        eta_eng = 0.3
    return m_oe, m_mto, m_f, pl_increase, LD_crs, eta_eng

def max_payload_mass(m_pldes, pl_increase, m_mto, m_oe):
    m_plmax = min(m_pldes * pl_increase, m_mto - m_oe)
    return m_plmax

def fuelmass_maxpl(m_mto, m_oe, m_plmax):
    m_fB = max(0, m_mto - m_oe - m_plmax)
    return m_fB

def R_cruise(eta_eng, LD_crs, m_oe, m_plmax, m_fuel):
    r_tot =  eta_eng * eta_p * (LD_crs) * (e_f/g) * np.log((m_oe + m_plmax + m_fuel)/(m_oe + m_plmax))
    return r_tot

def R_tot(R_nom, R_cruise, LD):
    R_lost = (1 / 0.7) * (LD) * (h_cr + ((V_cr ** 2) / (2 * g)))
    R_eq = ((R_nom + R_lost) * (1 + f_con)) + (1.2 * R_div) + (t_E * V_cr)
    R_aux = R_eq - R_nom
    R = R_cruise - R_aux
    return R

def plotting(ranges, plmasses, title, colour):
    # plotting the points
    plt.plot(ranges, plmasses, color=colour, linewidth=3,
             marker='o', markerfacecolor=colour, markersize=5)

    plt.axhline(y=plmasses[2], color='grey', linestyle='--')
    plt.annotate('Design payload', xy=(1700, plmasses[2] + 50))
    plt.axhline(y=plmasses[1], color='grey', linestyle='--')
    plt.annotate('Maximum payload', xy=(1700, plmasses[1] + 50))

    plt.axvline(x=ranges[2], color='grey', linestyle='--')
    plt.annotate('Range @ Design payload', xy=(ranges[2] - 60, 100), rotation='vertical')
    plt.axvline(x=ranges[1], color='grey', linestyle='--')
    plt.annotate('Range @ Maximum payload', xy=(ranges[1] - 60, 100), rotation='vertical')

    plt.xlim(0,2300)
    plt.ylim(0,7050)
    n = ['A', 'B', 'C', 'D']
    for i, txt in enumerate(n):
        plt.annotate(txt, (ranges[i], plmasses[i]))
    # naming the x axis
    plt.xlabel('Range [nmi]')
    # naming the y axis
    plt.ylabel('Payload [kg]')
    # giving a title to my graph
    plt.title(title)
    # function to show the plot
    plt.show()


# -------------------- LH2 FUEL CELL --------------------

m_oeFC, m_mtoFC, m_fFC, pl_increaseFC, LD_crsFC, eta_engFC = configuration_values(1)
m_plmaxFC = max_payload_mass(m_pldes, pl_increaseFC, m_mtoFC, m_oeFC)
m_fB_FC = fuelmass_maxpl(m_mtoFC, m_oeFC, m_plmaxFC)

# Point B
r_bFC = R_cruise(eta_engFC, LD_crsFC, m_oeFC,m_plmaxFC, m_fB_FC)
R_B_FC = R_tot(r_bFC, r_bFC, LD_crsFC)/1852
# Point C
r_cFC = R_cruise(eta_engFC, LD_crsFC, m_oeFC,m_pldes, m_fFC)
R_C_FC = R_tot(r_cFC, r_cFC, LD_crsFC)/1852
# Point D
r_dFC = R_cruise(eta_engFC, LD_crsFC, m_oeFC,0, m_fFC)
R_D_FC = R_tot(r_dFC, r_dFC, LD_crsFC)/1852

rangesFC = [0, R_B_FC, R_C_FC, R_D_FC]
plmassesFC = [m_plmaxFC, m_plmaxFC, m_pldes, 0]
print('LH2 Fuel cell', rangesFC, plmassesFC)

# -------------------- LH2 COMBUSTION --------------------

m_oeC, m_mtoC, m_fC, pl_increaseC, LD_crsC, eta_engC = configuration_values(2)
m_plmaxC = max_payload_mass(m_pldes, pl_increaseC, m_mtoC, m_oeC)
m_fB_C = fuelmass_maxpl(m_mtoC, m_oeC, m_plmaxC)

# Point B
r_bC = R_cruise(eta_engC, LD_crsC, m_oeC,m_plmaxC, m_fB_C)
R_B_C = R_tot(r_bC, r_bC, LD_crsC)/1852
# Point C
r_cC = R_cruise(eta_engC, LD_crsC, m_oeC,m_pldes, m_fC)
R_C_C = R_tot(r_cC, r_cC, LD_crsC)/1852
# Point D
r_dC = R_cruise(eta_engC, LD_crsC, m_oeC,0, m_fC)
R_D_C = R_tot(r_dC, r_dC, LD_crsC)/1852

rangesC = [0, R_B_C, R_C_C, R_D_C]
plmassesC = [m_plmaxC, m_plmaxC, m_pldes, 0]
print('LH2 Combustion', rangesC, plmassesC)
plotting(rangesFC, plmassesFC, 'Payload-Range diagram for LH2 fuel cell', 'red')
plotting(rangesC, plmassesC, 'Payload-Range diagram for LH2 Combustion', 'green')
