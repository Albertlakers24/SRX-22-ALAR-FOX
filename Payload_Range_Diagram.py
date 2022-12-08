#100% fuel calculations
import numpy as np
from matplotlib import pyplot as plt

n_pk =  0.85
n_ph = 0.85
PSFC = 0.48*(0.45/(745*3600))
m_pldes = 100*50
L_D = 16                                                      #Lift over drag
e_fk = 42.9 * 10**6
e_fh = 142 * 10**6
n_engk = (1/e_fk)*(1/PSFC)
n_engh = (1/e_fh)*(1/PSFC)
g = 9.81
V_cr = 143.89
h_cr = 11000
LD_crs = 16.7
#
m_oe =  12400                                                    #given                                                    #
#m_mto  =                                                    #given
# #Point A
ranges1 = [0]
plmasses1 = [m_pldes]
#Point B (Range at Max Payload)

R_b = n_eng * n_p * (L/D) * (e_f /g) * np.log((m_oe + m_plmax + m_f)/(m_oe + m_pl))
#
m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))
#
ranges = np.append(ranges, [R_b])
plmasses = np.append(plmasses, [m_pl])

print(ranges)

#Point C (Design Range)
f_con = 5/100
R_div = 200000
t_E = 45 * 60
R_lost1 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
print(R_lost1)
R_nom1 = 1852000
R_eq1 = ((R_nom1 + R_lost1)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
print(R_eq1)
R_aux1 = R_eq1 - R_nom1
print(R_aux1)
m_mto = 21500
m_f = m_mto * (1 - np.exp(-R_nom1 / (n_engh * n_ph * (e_fh /g) * (L_D))))
R_c1 = (n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_f)/(m_oe + m_pldes))) - R_aux1
print("Rc", R_c1)
ranges1 = np.append(ranges1, [R_c1])
plmasses1 = np.append(plmasses1, [m_pldes])

#Point D (Ferry Range)

R_d1 = (n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe  + m_f)/(m_oe ))) - R_aux1
print(R_d1)

ranges1 = np.append(ranges1, [R_d1])
plmasses1 = np.append(plmasses1, [0])


#Constructing the actual plot

# plotting the points
plt.plot(ranges1, plmasses1, color='green', linewidth=2,
         marker='o', markerfacecolor='green', markersize=5)


# naming the x axis
plt.xlabel('Range (m)')
# naming the y axis
plt.ylabel('Payload (Kg)')

# giving a title to my graph
plt.title('Payload-Range diagram for 100% combustion')

# function to show the plot
plt.show()

#------------

#Hydrogen- fuel case


# #Point A
ranges2 = [0]
plmasses2 = [m_pldes]
#Point B (Range at Max Payload)
# R_b = n_eng * n_p * (L/D) * (e_f /g) * np.log((m_oe + m_plmax + m_f)/(m_oe + m_pl))
#
# m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))
#

#m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))

ranges = np.append(ranges, [R_b])
plmasses = np.append(plmasses, [m_pl])

print(ranges)

#Point C (Design Range)
ratio_h = np.linspace(0.0, 1.0, num=21)
ratio_k = 1 - ratio_h
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
m_fh = ((m_mto * (1 - np.exp(-R_nom2 / (n_engh * n_ph * (e_fh /g) * (L_D)))))*0.3)
m_fk =  ((m_mto * (1 - np.exp(-R_nom2 / (n_engk * n_pk * (e_fk /g) * (L_D)))))*0.7)
R_c2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_fh)/(m_oe + m_pldes))) +  (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe + m_pldes + m_fk)/(m_oe + m_pldes))))   - R_aux1
print("Rc", R_c2)
ranges2 = np.append(ranges2, [R_c2])
plmasses2 = np.append(plmasses2, [m_pldes])

#Point D (Ferry Range)

R_d2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe  + m_fh)/(m_oe )))+ (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe  + m_fk)/(m_oe )))) - R_aux1
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

a = np.linspace(0.0, 1.0, num=21)
print(a)
b = 1-a

print(b)
