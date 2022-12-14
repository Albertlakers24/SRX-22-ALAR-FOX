import numpy as np
import matplotlib.pyplot as plt

g = 9.80665
prop_type = 1       # 1 = Parallel Series-Hybrid, 2 = Series

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time
m_pldes = 5440   #Design payload [kg]
LD_crs= 17

# Energy
e_f = 42*10**6
e_bat = 2.7*10**6
E_tot = 32604 * 10**6

# Efficiencies
eta_em = 0.95
eta_gt = 0.39
eta_p = 0.85
eta_gen = 0.97
eta_tp = 0.45
eta_bat = 0.934

# Degree of hybridisation
DoH_TO= 0.55
DoH_cruise = 0.2
t_TO = 381 + 1060 + 1981
t_cruise = 15624
t_descent = 720 + 541
phi = ((DoH_TO * t_TO) + (DoH_cruise*t_cruise))/(t_descent+t_TO+t_cruise)
print(phi)

phiB = 0.45

def configuration_values(prop_type):
    if prop_type == 1:  # Parallel Series
        m_mto = 24030
        m_oe = 12900
        m_f = 1720
        m_bat = 3970
        pl_increase = 1 #1.47
        eta_1 = eta_tp
        eta_2 = eta_em
        eta_3 = eta_p
        E_0tot = m_bat*e_bat + m_f*e_f

    if prop_type == 2:  # Series
        m_mto = 25290
        m_oe = 13700
        m_f = 2210
        m_bat = 3940
        pl_increase = 1 #1.54
        eta_1 =eta_gt*eta_gen
        eta_2 = eta_bat
        eta_3 = eta_em*eta_p
        E_0tot = m_bat * e_bat + m_f * e_f

    return m_mto, m_oe, m_f, m_bat, pl_increase, eta_1, eta_2, eta_3, E_0tot

def max_payload_mass(m_pldes, pl_increase, m_mto, m_oe, m_bat):
    m_plmax = min(m_pldes * pl_increase, m_mto - m_oe - m_bat)
    return m_plmax

def R_cruise(eta_1, eta_2, eta_3, LD, phi, m_oe, E_0tot, m_pl):
    r_tot = eta_3*(e_f/g) * (LD) *(eta_1 + eta_2 *(phi/(1-phi))) * np.log((m_oe + m_pl + ((g/e_bat)*E_0tot*(phi + ((e_bat/e_f)*(1 - phi)))))/(m_oe+m_pl+((g/e_bat)*phi*E_tot)))
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
    plt.ylim(0,8700)
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

# ----------------------------- PARALLEL SERIES -----------------------------
m_mtoPS, m_oePS, m_fPS, m_batPS, pl_increasePS, eta_1PS, eta_2PS, eta_3PS, E_0totPS = configuration_values(1)
m_plmaxPS = max_payload_mass(m_pldes, pl_increasePS, m_mtoPS, m_oePS, m_batPS)

# Point B
R_b = R_cruise(eta_1PS, eta_2PS, eta_3PS, LD_crs, phiB, m_oePS, E_0totPS, m_plmaxPS)  # (eta_1, eta_2, eta_3, LD, phi, m_oe, E_tot, m_pl)
RB = R_tot(R_b, R_b, LD_crs)/1852       # (R_nom, R_cruise, LD)
# Point C
R_c = R_cruise(eta_1PS, eta_2PS, eta_3PS, LD_crs, phi, m_oePS, E_0totPS, m_pldes)  # (eta_1, eta_2, eta_3, LD, phi, m_oe, E_tot, m_pl)
RC = R_tot(R_c, R_c, LD_crs)/1852       # (R_nom, R_cruise, LD)
# Point D
R_d = R_cruise(eta_1PS, eta_2PS, eta_3PS, LD_crs, phi, m_oePS, E_0totPS, 0)  # (eta_1, eta_2, eta_3, LD, phi, m_oe, E_tot, m_pl)
RD = R_tot(R_d, R_d, LD_crs)/1852       # (R_nom, R_cruise, LD)

rangesPS = [0, RB, RC, RD]
plmassesPS = [m_plmaxPS, m_plmaxPS,m_pldes, 0]
print('Hybrid Electric Parallel Series', rangesPS, plmassesPS)

# ----------------------------- SERIES CONFIGURATION -----------------------------
m_mtoS, m_oeS, m_fS, m_batS, pl_increaseS, eta_1S, eta_2S, eta_3S, E_0totS = configuration_values(2)
m_plmaxS = max_payload_mass(m_pldes, pl_increaseS, m_mtoS, m_oeS, m_batS)

# Point B
R_b = R_cruise(eta_1S, eta_2S, eta_3S, LD_crs, phiB, m_oeS, E_0totS, m_plmaxS)  # (eta_1, eta_2, eta_3, LD, phi, m_oe, E_tot, m_pl)
RB = R_tot(R_b, R_b, LD_crs)/1852       # (R_nom, R_cruise, LD)
# Point C
R_c = R_cruise(eta_1S, eta_2S, eta_3S, LD_crs, phi, m_oeS, E_0totS, m_pldes)  # (eta_1, eta_2, eta_3, LD, phi, m_oe, E_tot, m_pl)
RC = R_tot(R_c, R_c, LD_crs)/1852       # (R_nom, R_cruise, LD)
# Point D
R_d = R_cruise(eta_1S, eta_2S, eta_3S, LD_crs, phi, m_oeS, E_0totS, 0)  # (eta_1, eta_2, eta_3, LD, phi, m_oe, E_tot, m_pl)
RD = R_tot(R_d, R_d, LD_crs)/1852       # (R_nom, R_cruise, LD)

rangesS = [0, RB, RC, RD]
plmassesS = [m_plmaxS, m_plmaxS,m_pldes, 0]
print('Hybrid Electric Series', rangesS, plmassesS)

plotting(rangesS, plmassesS, 'Payload range diagram Hybrid Electric Series', 'blue')
plotting(rangesPS, plmassesPS, 'Payload range diagram Hybrid Electric Parallel Series', 'orange')

'''
WHY THIS METHOD SUCKS:
The results of the New VOS method are incoherent relative to similar aircraft, due to the fact that it uses a
constant average DoH throughout the entire flight. As such, this method merely proposes an incorrect estimation, 
as the average DoH will result in power values that are far in excess of the required during certain flight phases,
while also providing power values that are below the required for other flight phases. As a result, this method is
incorrect for a significant portion of the flight, and therefore the NEW method is superior in this regard. All
calculations will, from this point forward, be performed with the NEW method.'''
