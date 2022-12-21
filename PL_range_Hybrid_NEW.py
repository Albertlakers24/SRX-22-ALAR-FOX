import numpy as np
from matplotlib import pyplot as plt

# constants
g = 9.80665
e_bat       = 1.656  * 10**6
e_f         = 42 * 10**6

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time

# Common values for both configurations
m_pldes = 5443                    #Design payload [kg]
E_tot =25365 * 10**6             # Total propulsive energy (in J)
LD_crs     = 22

# Efficiencies
eta_p = 0.85                   # Propulsive efficiency (overall)
eta_i = 0.99                        # Inverter efficiency
#n_eng_em= 0.934*0.995*0.95 # 0.934*0.99*0.995*0.95*0.85 #Enine efficiency (electric motor)
eta_m = 0.95
c_b = e_bat                         # Battery specific energy (J/kg)

def configuration_values(prop_type):
    if prop_type == 1:  # Parallel Series
        m_mto = 33346
        m_oe = 19446
        m_f = 2159
        m_bat = 6298
        pl_increase = 1.2
        c_p = 1/(43*10**6*0.45)
        eta_g = 0.45  # Generator efficiency

    if prop_type == 2:  # Series
        m_mto = 38532
        m_oe = 21785
        m_f = 2778
        m_bat = 8525
        pl_increase = 1.25
        c_p = 1/(43*10**6*0.39)
        eta_g = 0.97*0.99*0.99*0.39  # Generator efficiency

    return m_mto, m_oe, m_f, m_bat, pl_increase, c_p, eta_g

def max_payload_mass(m_pldes, pl_increase, m_mto, m_oe, m_bat):
    m_plmax = min(m_pldes * pl_increase, m_mto - m_oe - m_bat)
    return m_plmax

def fuelmass_maxpl(m_mto, m_oe, m_bat, m_plmax):
    m_fB = max(0, m_mto - m_oe - m_bat - m_plmax)
    return m_fB

def f_ratio(m_mto, m_oe):
    f = m_oe/m_mto
    return f

def mass_fraction(m_mto,m_bat, m_fuel, f):
    phi_bat = m_bat / m_mto / (1 - f)
    psi_fuel = m_fuel / m_mto / (1 - f)
    return phi_bat, psi_fuel

def R_cruise(phi, psi, f, c_p, eta_g):
    r_tot = LD_crs * eta_i * eta_m * eta_p / g * ((phi * c_b * (1 - f)) / (1 - psi * (1 - f)) + eta_g / c_p * np.log(1 / (1 - psi * (1 - f))))
    return r_tot

def R_tot(R_nom, R_cruise, LD_crs):
    R_lost = (1 / 0.7) * (LD_crs) * (h_cr + ((V_cr ** 2) / (2 * g)))
    R_eq = ((R_nom + R_lost) * (1 + f_con)) + (1.2 * R_div) + (t_E * V_cr)
    R_aux = R_eq - R_nom
    R = R_cruise + R_aux
    return R_eq

def plotting(ranges, plmasses, title, colour):
    # plotting the points
    plt.plot(ranges, plmasses, color=colour, linewidth=3,
             marker='o', markerfacecolor=colour, markersize=5)

    plt.axhline(y=plmasses[2], color='grey', linestyle='--')
    plt.annotate('Design payload', xy=(1700, plmasses[2] + 50))
    plt.axhline(y=plmasses[1], color='grey', linestyle='--')
    plt.annotate('Maximum payload', xy=(1700, plmasses[1] + 50))

    plt.axvline(x=ranges[2], color='grey', linestyle='--')
    plt.annotate('Range @ Design payload', xy=(ranges[2] - 60, 200), rotation='vertical')
    plt.axvline(x=ranges[1], color='grey', linestyle='--')
    plt.annotate('Range @ Maximum payload', xy=(ranges[1] - 60, 200), rotation='vertical')

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


# ----------------------------- PARALLEL SERIES CONFIGURATION -----------------------------
m_mtoPS, m_oePS, m_fPS, m_batPS, pl_increasePS, c_pPS, eta_g_PS = configuration_values(1)
m_plmaxPS = max_payload_mass(m_pldes, pl_increasePS, m_mtoPS, m_oePS, m_batPS)
# Point B
fPS = f_ratio(m_mtoPS, m_oePS)
m_fB_PS = fuelmass_maxpl(m_mtoPS, m_oePS, m_batPS, m_plmaxPS)
phi_B_PS, psi_B_PS = mass_fraction(m_mtoPS, m_batPS, m_fB_PS, fPS)
R_b_PS = R_cruise(phi_B_PS, psi_B_PS, fPS, c_pPS,eta_g_PS)
RB_PS = R_tot(R_b_PS, R_b_PS, LD_crs)/1852
# Point C
phi_C_PS, psi_C_PS = mass_fraction(m_mtoPS, m_batPS, m_fPS, fPS)
R_c_PS = R_cruise(phi_C_PS, psi_C_PS,fPS, c_pPS,eta_g_PS)
RC_PS = R_tot(R_c_PS, R_c_PS, LD_crs)/1852
# Point D
m_mtoD_PS = m_oePS + m_batPS + m_fPS
f_D_PS = f_ratio(m_mtoD_PS, m_oePS)
phi_D_PS, psi_D_PS = mass_fraction(m_mtoD_PS, m_batPS, m_fPS, f_D_PS)
R_d_PS = R_cruise(phi_D_PS, psi_D_PS,f_D_PS, c_pPS,eta_g_PS)
RD_PS = R_tot(R_d_PS, R_d_PS, LD_crs)/1852

rangesPS = [0, RB_PS, RC_PS, RD_PS]
plmassesPS = [m_plmaxPS, m_plmaxPS,m_pldes, 0]
print('Hybrid Electric Parallel Series', rangesPS, plmassesPS)

# -----------------------------  SERIES CONFIGURATION -----------------------------

m_mtoS, m_oeS, m_fS, m_batS, pl_increaseS, c_pS, eta_g_S = configuration_values(2)
m_plmaxS = max_payload_mass(m_pldes, pl_increaseS, m_mtoS, m_oeS, m_batS)
# Point B
fS = f_ratio(m_mtoS, m_oeS)
m_fB_S = fuelmass_maxpl(m_mtoS, m_oeS, m_batS, m_plmaxS)
phi_B_S, psi_B_S = mass_fraction(m_mtoS, m_batS, m_fB_S, fS)
R_b_S = R_cruise(phi_B_PS, psi_B_PS, fS, c_pS,eta_g_S)
RB_S = R_tot(R_b_S, R_b_S, LD_crs)/1852
# Point C
phi_C_S, psi_C_S = mass_fraction(m_mtoS, m_batS, m_fS, fS)
R_c_S = R_cruise(phi_C_S, psi_C_S, fS, c_pS,eta_g_S)
RC_S = R_tot(R_c_S, R_c_S, LD_crs)/1852
# Point D
m_mtoD_S = m_oeS + m_batS + m_fS
f_D_S = f_ratio(m_mtoD_S, m_oeS)
phi_D_S, psi_D_S = mass_fraction(m_mtoD_S, m_batS, m_fS, f_D_S)
R_d_S = R_cruise(phi_D_S, psi_D_S,fS, c_pS,eta_g_S)
RD_S = R_tot(R_d_S, R_d_S, LD_crs)/1852

rangesS = [0, RB_S, RC_S, RD_S]
plmassesS = [m_plmaxS, m_plmaxS,m_pldes, 0]
print('Hybrid Electric Series',rangesS, plmassesS)


plotting(rangesS, plmassesS, 'Payload range diagram Electric Hybrid Series', 'blue')
plotting(rangesPS, plmassesPS, 'Payload range diagram Hybrid Electric Parallel Series', 'orange')
