import numpy as np
from matplotlib import pyplot as plt

g = 9.80665
prop_type = 1       # 1 = Parallel Series-Hybrid, 2 = Series

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time

m_pldes = 5443                    #Design payload [kg]
m_plmax = m_pldes*1.2
m_bat       = 3898                      #battery mass [kg]
E_tot = 31226 * 10**6                 # Total propulsive energy (in J)

#Aircraft mass characteristics -- --> to be updated!!
#PARALLEL SERIES
# m_oe    =  12883                    #Operational empty weight [kg]
# m_mto   = 25857                    # Max take-off
# phi = (0.35+0.02)/2                               # Rate of hybridization (Point C - Design point)
# n_eng_tp = 0.45          #Engine efficiency (thermodynamic, turboprop)
# m_fuel = 3633
# phiB = (0.45+0.02)/2

# SERIES
n_eng_tp = 0.39*0.97*0.99*0.99
m_oe    =  12786                   # Operational empty weight [kg]
m_mto   = 24378                    # Max take-off
phi = (0.45+ 0.04)/2               # Rate of hybridization (Point C - Design point)
m_fuel = 2211
phiB = (0.65+0.04)/2

#Propulsion characteristics

# ----- BATTERY ------
e_bat       = 2.7 * 10**6
n_eng_em    = 0.934*0.85*0.99*0.995*0.95                      #Enine efficiency (electric motor)
n_p_em      = 0.85                      #Propulsive efficiency (electric motor)
#P_em = 2000*10**3
# ------- TURBOPROP
#PSFC        = 0.48*(0.45/(745*3600))    #Specific fuel consumption
e_f         = 43 * 10**6


n_p_tp      = 0.85                      #Propulsive efficiency (turboprop)
#P_tp = 8000*10**3
# ----- HYBRID THINGS
n_p = 0.85                              # Propulsive efficiency (overall)

#---- PHI Calculations for Point B (Point C included for sanity check)
m_f = m_mto - m_oe - m_pldes - m_bat
# print(m_fuel, m_f)
# phi1 =(m_f*e_f/E_tot/g)
# print(phi, phi1)
m_fB = m_mto - m_oe - m_plmax - m_bat       # Mass of fuel at max payload - battery mass constant - trading fuel for payload
# phiB =  (m_fB*e_f/E_tot/g)               # Rate of hybridization (Point B - Max payload point)
# print(m_fB, phiB)


#aerodynamic characteristics
LD     = 17
LD_crs = 16.7              #Lift over drag (cruise)

#tbs
# m_f     = 0.8 * (m_mto * (1 - np.exp(-R_nom / (n_eng_tp * n_p_tp * (e_f /g) * (LD_crs)))))
# print(m_f)


#Design Point A

#Design Point B
R_fB = (n_eng_tp * n_p_tp * (LD_crs) * (e_f /g) *  np.log((m_oe + m_plmax + m_fB) / (m_oe + m_plmax)))        #Brequet range equation (fuel)
R_bB = (n_eng_em * n_p_em * (LD_crs) * (e_bat /g) * ((m_bat / (m_oe + m_plmax + m_bat)))     )                 #Brequet range equation (battery)

RB = n_p*(e_f/g) * (LD_crs) *(n_eng_tp + n_eng_em *(phiB/(1+phiB))) * np.log((m_oe + m_plmax + ((g/e_bat)*E_tot*(phiB + ((e_bat/e_f)*(1 - phiB)))))/(m_oe+m_plmax+((g/e_bat)*phiB*E_tot)))

R_lostB3 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eqB3 = ((R_fB + R_bB+ R_lostB3)*(1+f_con))  + (1.2*R_div) + (t_E * V_cr)
R_auxB3 = R_eqB3 - R_fB - R_bB
R_b = R_fB + R_bB - R_auxB3
R_B2 = RB - R_auxB3


#point C (design range)
#range equations
m_f = m_mto - m_oe - m_bat - m_pldes
R_lost          = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eq            = ((R_nom + R_lost)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)                                         #Equivalent range
R_aux           = R_eq - R_nom                                                                                      #Auxilary range
R_f             = (n_eng_tp * n_p_tp * (LD_crs) * (e_f /g) *  np.log((m_oe + m_pldes + m_f) / (m_oe + m_pldes)))        #Brequet range equation (fuel)
R_b             = (n_eng_em * n_p_em * (LD_crs) * (e_bat /g) * ((m_bat / (m_oe + m_pldes + m_bat)))     )                 #Brequet range equation (battery)


R_c             = R_f + R_b - R_aux

RC = n_p*(e_f/g) * (LD_crs) *(n_eng_tp + n_eng_em *(phi/(1+phi))) * np.log((m_oe + m_pldes + ((g/e_bat)*E_tot*(phi + ((e_bat/e_f)*(1 - phi)))))/(m_oe+m_pldes+((g/e_bat)*phi*E_tot)))
R_C2 = RC - R_aux


#Point D (Ferry Range)
R_f_ferry       = (n_eng_tp * n_p_tp * (LD_crs) * (e_f /g) *  np.log((m_oe + m_f) / (m_oe)))                            #Brequet range equation (fuel)
R_b_ferry       = (n_eng_em * n_p_em * (LD_crs) * (e_bat /g) * ((m_bat / (m_oe + m_bat)))     )                           #Brequet range equation (battery)
R_d = R_f_ferry + R_b_ferry - R_aux

RD = n_p*(e_f/g) * (LD_crs) *(n_eng_tp + n_eng_em *(phi/(1+phi))) * np.log((m_oe + ((g/e_bat)*E_tot*(phi + ((e_bat/e_f)*(1 - phi)))))/(m_oe+((g/e_bat)*phi*E_tot)))
R_D2 = RD - R_aux


ranges = np.array([0, R_B2/1852, R_C2/1852, R_D2/1852])
plmasses = np.array([m_plmax, m_plmax, m_pldes,0])

print(ranges, plmasses)

#Constructing the actual plot

# plotting the points
plt.plot(ranges, plmasses, color='orange', linewidth=3,
         marker='o', markerfacecolor='orange', markersize=5)
n = ['A','B','C','D']
for i, txt in enumerate(n):
    plt.annotate(txt, (ranges[i], plmasses[i]))

# naming the x axis
plt.xlabel('Range [nmi]')
# naming the y axis
plt.ylabel('Payload [kg]')

# giving a title to my graph
plt.title('Payload range diagram Electric-Hybrid Series')

# function to show the plot
plt.show()



