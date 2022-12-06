import numpy as np
from matplotlib import pyplot as plt

n_pk =  0.85                                         #Propulsive efficiency twin turboprop
n_ph = 0.85                                          #Propulsive efficiency twin turboprop
PSFC = 0.48*(0.45/(745*3600))                        #
m_pldes = 100*50                                     # Payload mass designed for
L_D = 16                                             #Lift to drag ratio
e_fk = 42.9 * 10**6
e_fh = 142 * 10**6
n_engk = (1/e_fk)*(1/PSFC)
n_engh = (1/e_fh)*(1/PSFC)
g = 9.81
V_cr = 143.89                                         # Cruise speed in m/s
h_cr = 11000                                          # Cruise altitude in m
LD_crs = 16.7                                         # Lift to drag ratio during cruise
#
m_oe =  12400                                         #Operating empty mass

# #Point A
ranges2 = [0]
plmasses2 = [m_pldes]
#Point B (Range at Max Payload)
# m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))
#
# R_b = n_eng * n_p * (L/D) * (e_f /g) * np.log((m_oe + m_plmax + m_f)/(m_oe + m_pl))
#
# m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))
#
# ranges = np.append(ranges, [R_b])
# plmasses = np.append(plmasses, [m_pl])

#print(ranges)

#Point C (Design Range)
ratio_h = 0.3
ratio_k = 0.7
f_con = 5/100
R_div = 200000
t_E = 45 * 60
R_lost2 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
print(R_lost2)
R_nom2 = 1852000
R_eq2 = ((R_nom2 + R_lost2)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
print(R_eq2)
R_aux2 = R_eq2 - R_nom2
print(R_aux2)
m_mto = 21500
m_fh = ((m_mto * (1 - np.exp(-R_nom2 / (n_engh * n_ph * (e_fh /g) * (L_D)))))*ratio_h)
m_fk =  ((m_mto * (1 - np.exp(-R_nom2 / (n_engk * n_pk * (e_fk /g) * (L_D)))))*ratio_k)
R_c2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_fh)/(m_oe + m_pldes))) +  (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe + m_pldes + m_fk)/(m_oe + m_pldes))))   - R_aux2
print("Rc", R_c2)
ranges2 = np.append(ranges2, [R_c2])
plmasses2 = np.append(plmasses2, [m_pldes])

#Point D (Ferry Range)

R_d2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe  + m_fh)/(m_oe )))+ (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe  + m_fk)/(m_oe )))) - R_aux2
print(R_d2)

ranges2 = np.append(ranges2, [R_d2])
plmasses2 = np.append(plmasses2, [0])

#Constructing the actual plot

# plotting the points
plt.plot(ranges2, plmasses2, color='red', linewidth=2,
         marker='o', markerfacecolor='red', markersize=5)


# naming the x axis
plt.xlabel('Range (m)')
# naming the y axis
plt.ylabel('Payload (Kg)')

# giving a title to my graph
plt.title('Payload-Range diagram for LH-Kerosene hybrid')

# function to show the plot
plt.show()

print()