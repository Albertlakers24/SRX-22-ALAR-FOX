#100% fuel calculations
import numpy as np
from matplotlib import pyplot as plt

n_p =  0.85                                                      #Propulsive efficiency
PSFC = 0.48*(0.45/(745*3600))
m_pldes = 100*50
L_D = 16                                                      #Lift over drag
e_f = 42.8 * 10**6
n_eng = (1/e_f)*(1/PSFC)
g = 9.81
V_cr = 143.89
h_cr = 11000
LD_crs = 16.7
#
m_oe =  12400                                                    #given                                                    #
#m_mto  =                                                    #given
# #Point A
ranges = [0]
plmasses = [m_pldes]
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
f_con = 5/100
R_div = 200000
t_E = 45 * 60
R_lost = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
print(R_lost)
R_nom = 1852000
R_eq = ((R_nom + R_lost)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
print(R_eq)
R_aux = R_eq - R_nom
print(R_aux)
m_mto = 21500
m_f = m_mto * (1 - np.exp(-R_nom / (n_eng * n_p * (e_f /g) * (L_D))))
R_c = (n_eng * n_p * (L_D) * (e_f /g) * np.log((m_oe + m_pldes + m_f)/(m_oe + m_pldes))) - R_aux
print("Rc", R_c)
ranges = np.append(ranges, [R_c])
plmasses = np.append(plmasses, [m_pldes])

#Point D (Ferry Range)

R_d = R_c = (n_eng * n_p * (L_D) * (e_f /g) * np.log((m_oe  + m_f)/(m_oe ))) - R_aux
print(R_d)

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
plt.title('Some cool customizations!')

# function to show the plot
plt.show()

