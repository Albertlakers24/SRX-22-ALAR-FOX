import numpy as np
from matplotlib import pyplot as plt

prop_type = 2 # 1 for Parallel Series, 2 for Series

# constants
g = 9.80665
e_bat       = 2.7 * 10**6
e_f         = 42 * 10**6

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time

# Common values for both confugurations
m_pldes = 5443                    #Design payload [kg]
E_tot = 32604 * 10**6             # Total propulsive energy (in J)
LD     = 17
n_eng_em= 0.934*0.85*0.99*0.995*0.95   #Enine efficiency (electric motor)


eta_p = 0.85                        # Propulsive efficiency (overall)
eta_i = 0.99
eta_m = n_eng_em
eta_g = 0.97
c_b = e_bat
c_p  = 0.48*(0.45/(745*3600))  # Specific fuel consumption
c_p = 300/3600/1000/1000
print(c_p)
BSFC= 1/(43*10**6 * 0.45)      #1/(43*10**6 * 0.39 * 0.9 *0.99)   #Brake-specific fuel consumption (only 43*10^6 * 0.45 if parallel series)
print(BSFC)
if prop_type == 1:
    m_mto = 24100
    m_oe = 12883
    m_f = 1723
    m_bat = 3970
    m_plmax = m_pldes * 1.535

if prop_type == 2:
    m_mto = 25300
    m_oe = 13700
    m_f = 2211
    m_bat = 3940
    m_plmax = m_pldes * 1.54

m_fB = m_mto - m_oe - m_bat - m_plmax
m_mtoD = m_oe + m_bat + m_f
f = m_oe / m_mto
f_D = m_oe/m_mtoD
print(m_mto)
def mass_fraction(m_mto,m_bat, m_fuel, f):
    phi_bat = m_bat / m_mto / (1 - f)
    psi_fuel = m_fuel / m_mto / (1 - f)
    return phi_bat, psi_fuel

def W_pay(phi,psi):
	w_pay = (1-f)*m_mto*(1-phi-psi)
	return w_pay

def R_cruise(phi, psi):
    r_tot = LD * eta_i * eta_m * eta_p / g * ((phi * c_b * (1 - f)) / (1 - psi * (1 - f)) + eta_g / c_p * np.log(1 / (1 - psi * (1 - f))))
    return r_tot

def R_tot(R_nom, R_cruise):
    R_lost = (1 / 0.7) * (LD) * (h_cr + ((V_cr ** 2) / (2 * g)))
    R_eq = ((R_nom + R_lost) * (1 + f_con)) + (1.2 * R_div) + (t_E * V_cr)
    R_aux = R_eq - R_nom
    R = R_cruise + R_aux
    return R

def plotting(ranges, plmasses, title, colour):
    # plotting the points
    plt.plot(ranges, plmasses, color=colour, linewidth=3,
             marker='o', markerfacecolor=colour, markersize=5)

    plt.axhline(y=plmasses[2], color='grey', linestyle='--')
    plt.annotate('Design payload', xy=(1300, plmasses[2] + 50))
    plt.axhline(y=plmasses[1], color='grey', linestyle='--')
    plt.annotate('Maximum payload', xy=(1300, plmasses[1] + 50))

    plt.axvline(x=ranges[2], color='grey', linestyle='--')
    plt.annotate('Range @ Design payload', xy=(ranges[2] - 60, 100), rotation='vertical')
    plt.axvline(x=ranges[1], color='grey', linestyle='--')
    plt.annotate('Range @ Maximum payload', xy=(ranges[1] - 60, 100), rotation='vertical')

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


# Calculations ----> Point B
phi_B, psi_B = mass_fraction(m_mto, m_bat, m_fB, f)
R_b = R_cruise(phi_B, psi_B)
RB = R_tot(R_b, R_b)/1852
# Calculations -----> Point C
phi_C, psi_C = mass_fraction(m_mto, m_bat, m_f, f)
R_c = R_cruise(phi_C, psi_C)
RC = R_tot(R_nom, R_c)/1852
#Calculations -----> Point D
phi_D, psi_D = mass_fraction(m_mtoD, m_bat, m_f, f_D)
R_d = R_cruise(phi_D, psi_D)
RD = R_tot(R_d, R_d)/1852

ranges = [0, RB, RC, RD]
plmasses = [m_plmax, m_plmax,m_pldes, 0]
print(ranges, plmasses)
if prop_type == 1:
    plotting(ranges, plmasses, 'Payload range diagram Electric Hybrid Parallel Series', 'orange')
if prop_type == 2:
    plotting(ranges, plmasses, 'Payload range diagram Electric Hybrid Series', 'blue')

