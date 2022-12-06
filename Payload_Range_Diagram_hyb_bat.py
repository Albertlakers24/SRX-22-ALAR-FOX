import numpy as np
from matplotlib import pyplot as plt

g = 9.81

#mission characteristics
h_cr = 11000                        #Cruise altitude [m]
V_cr = 143.89                       #Cruise speed [m/s]
R_nom = 1852000                     #Nominal Range
f_con = 5/100                       #fuel trip fraction for contingency
R_div = 200000                      #Fuel needed for missed approach
t_E = 45 * 60                       #Endurance time

#Aircraft mass characteristics
m_oe    =  12400                    #Operational empty weight (reference for now) [kg]
m_pldes = 100*50                    #Design payload [kg]
m_mto   = 21500

#Propulsion characteristics
m_bat       = 5000                      #battery mass [kg]
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

ranges = [0]
plmasses = [m_pldes]



#range equations
R_lost          = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eq            = ((R_nom + R_lost)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)                                         #Equivalent range
R_aux           = R_eq - R_nom                                                                                      #Auxilary range
R_f             = (n_eng_tp * n_p_tp * (LD) * (e_f /g) *  np.log((m_oe + m_pldes + m_f) / (m_oe + m_pldes)))        #Brequet range equation (fuel)
R_f_ferry       = (n_eng_tp * n_p_tp * (LD) * (e_f /g) *  np.log((m_oe + m_f) / (m_oe)))                            #Brequet range equation (fuel)
R_b             = (n_eng_em * n_p_em * (LD) * (e_bat /g) * ((m_bat / (m_oe + m_pldes + m_bat)))     )                 #Brequet range equation (battery)
R_b_ferry       = (n_eng_em * n_p_em * (LD) * (e_bat /g) * ((m_bat / (m_oe + m_bat)))     )                           #Brequet range equation (battery)

R_c             = R_f + R_b - R_aux

#point C (design range)
ranges = np.append(ranges, [R_c])
plmasses = np.append(plmasses, [m_pldes])

#Point D (Ferry Range)

R_d = R_f_ferry + R_b_ferry - R_aux


ranges = np.append(ranges, [R_d])
plmasses = np.append(plmasses, [0])


#Constructing the actual plot

# plotting the points
plt.plot(ranges, plmasses, color='green', linewidth=3,
         marker='o', markerfacecolor='blue', markersize=12)


# naming the x axis
plt.xlabel('Range')
# naming the y axis
plt.ylabel('Payload')

# giving a title to my graph
plt.title('Payload range diagram hybrid battery/fuel')

# function to show the plot
plt.show()



