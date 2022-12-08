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

#Aircraft mass characteristics
m_oe    =  15454                    #Operational empty weight (reference for now) [kg]
m_pldes = 8124                    #Design payload [kg]
m_mto   = 30000
m_plmax = m_pldes*1.25

#Propulsion characteristics
m_bat       = 2600                      #battery mass [kg]
e_f         = 42.8 * 10**6
e_bat       = 2.7 * 10**6
PSFC        = 0.48*(0.45/(745*3600))    #Specific fuel consumption
n_eng_tp    = (1/e_f)*(1/PSFC)          #Engine efficiency (thermodynamic, turboprop)
n_eng_em    = 0.95                      #Enine efficiency (electric motor)
n_p_tp      = 0.85                      #Propulsive efficiency (turboprop)
n_p_em      = 0.95                      #Propulsive efficiency (electric motor)

#aerodynamic characteristics
LD     = 17
LD_crs = 16.7              #Lift over drag (cruise)

#tbs
m_f     = 0.8 * (m_mto * (1 - np.exp(-R_nom / (n_eng_tp * n_p_tp * (e_f /g) * (LD)))))


#Lists

#Design Point A
ranges3 = [0]
plmasses3 = [m_plmax]

#Design Point B
m_fB = m_mto - m_oe - m_plmax - m_bat
R_fB = (n_eng_tp * n_p_tp * (LD) * (e_f /g) *  np.log((m_oe + m_plmax + m_fB) / (m_oe + m_plmax)))        #Brequet range equation (fuel)
R_bB = (n_eng_em * n_p_em * (LD) * (e_bat /g) * ((m_bat / (m_oe + m_plmax + m_bat)))     )                 #Brequet range equation (battery)

R_lostB3 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eqB3 = ((R_fB + R_bB+ R_lostB3)*(1+f_con))  + (1.2*R_div) + (t_E * V_cr)
R_auxB3 = R_eqB3 - R_fB - R_bB
R_b = R_fB + R_bB - R_auxB3

ranges3 = np.append(ranges3, [R_b])
plmasses3 = np.append(plmasses3, [m_plmax])


#point C (design range)
#range equations
R_lost          = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eq            = ((R_nom + R_lost)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)                                         #Equivalent range
R_aux           = R_eq - R_nom                                                                                      #Auxilary range
R_f             = (n_eng_tp * n_p_tp * (LD) * (e_f /g) *  np.log((m_oe + m_pldes + m_f) / (m_oe + m_pldes)))        #Brequet range equation (fuel)
R_b             = (n_eng_em * n_p_em * (LD) * (e_bat /g) * ((m_bat / (m_oe + m_pldes + m_bat)))     )                 #Brequet range equation (battery)

R_c             = R_f + R_b - R_aux

ranges3 = np.append(ranges3, [R_c])
plmasses3 = np.append(plmasses3, [m_pldes])

#Point D (Ferry Range)
R_f_ferry       = (n_eng_tp * n_p_tp * (LD) * (e_f /g) *  np.log((m_oe + m_f) / (m_oe)))                            #Brequet range equation (fuel)
R_b_ferry       = (n_eng_em * n_p_em * (LD) * (e_bat /g) * ((m_bat / (m_oe + m_bat)))     )                           #Brequet range equation (battery)
R_d = R_f_ferry + R_b_ferry - R_aux


ranges3 = np.append(ranges3, [R_d])
plmasses3 = np.append(plmasses3, [0])


#Constructing the actual plot

# plotting the points
plt.plot(ranges3, plmasses3, color='green', linewidth=3,
         marker='o', markerfacecolor='blue', markersize=12)


# naming the x axis
plt.xlabel('Range')
# naming the y axis
plt.ylabel('Payload')

# giving a title to my graph
plt.title('Payload range diagram hybrid battery/fuel')

# function to show the plot
plt.show()



