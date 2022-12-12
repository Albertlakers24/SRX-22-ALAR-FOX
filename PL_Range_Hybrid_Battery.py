import numpy as np
from matplotlib import pyplot as plt

g = 9.80665

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time

#Aircraft mass characteristics -- --> to be updated!!
m_oe    =  15454                    #Operational empty weight [kg]
m_pldes = 5730                    #Design payload [kg]
m_mto   = 30000                    # Max take-off
m_plmax = m_pldes*1.3

#Propulsion characteristics

# ----- BATTERY ------
m_bat       = 2600                      #battery mass [kg]
e_bat       = 2.7 * 10**6
n_eng_em    = 0.95                      #Enine efficiency (electric motor)
n_p_em      = 0.95                      #Propulsive efficiency (electric motor)
P_em = 2000*10**3
# ------- TURBOPROP
PSFC        = 0.48*(0.45/(745*3600))    #Specific fuel consumption
e_f         = 42.8 * 10**6
n_eng_tp    = (1/e_f)*(1/PSFC)          #Engine efficiency (thermodynamic, turboprop)
n_p_tp      = 0.85                      #Propulsive efficiency (turboprop)
P_tp = 8000*10**3
# ----- HYBRID THINGS
E_tot = 35485.5 * 10**6                 # Total propulsive energy (in J)
n_p = 0.85                              # Propulsive efficiency (overall)
phi = 0.3                               # Rate of hybridization (Point C - Design point)

#---- PHI Calculations for Point B (Point C included for sanity check)
m_f = m_mto - m_oe - m_pldes - m_bat
phi1 = 1 - (m_f*e_f/E_tot/g)
print(phi, phi1)
m_fB = m_mto - m_oe - m_plmax - m_bat       # Mass of fuel at max payload - battery mass constant - trading fuel for payload
phiB = 1 - (m_fB*e_f/E_tot/g)               # Rate of hybridization (Point B - Max payload point)
print(phiB)

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
plt.plot(ranges, plmasses, color='blue', linewidth=3,
         marker='o', markerfacecolor='blue', markersize=5)
n = ['A','B','C','D']
for i, txt in enumerate(n):
    plt.annotate(txt, (ranges[i], plmasses[i]))

# naming the x axis
plt.xlabel('Range [nmi]')
# naming the y axis
plt.ylabel('Payload [kg]')

# giving a title to my graph
plt.title('Payload range diagram hybrid battery/fuel')

# function to show the plot
plt.show()



