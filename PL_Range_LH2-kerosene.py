import numpy as np
from matplotlib import pyplot as plt

prop_system = 2  # put 1 for LH2_Kerosene fuel cell, 2 for LH2 fuel cell

g = 9.80665

# Masses Aircraft --> Payload - constant for all configurations

m_pldes = 5730  # Payload mass designed for
m_plmax = m_pldes*1.1

# Masses Aircraft --> Max, empty, fuel - per propulsion system
if prop_system == 1:   # LH2_Kerosene fuel cell
    m_oe = 14643  # Operating empty mass
    m_mto = 21207
    m_fk = 390 # Kerosene mass in kg
    m_fh = 445 # LH2 mass in kg
if prop_system == 2:  # LH2 fuel cell
    m_oe = 14175  # Operating empty mass
    m_mto = 20434
    m_fh = 530 # LH2 mass in kg
    m_fk = 0   # Kerosene mass in kg

# Propulsion charecteristics
n_pk =  0.85                            #Propulsive efficiency twin turboprop
n_ph = 0.85                             #Propulsive efficiency twin turboprop
e_fk = 42.9 * 10**6                     # Specific Energy - kerosene
e_fh = 142 * 10**6                      # specific energy - hydrogen

if prop_system == 1:    # LH2_Kerosene fuel cell
    n_engk = 0.39*0.9*0.99                  # Efficiency engine - kerosene in fuel cell system
    n_engh = 0.6*0.97*(0.995**2)*0.85*0.95       # Efficiency engine - hydrogen fuel cell w kerosene
if prop_system == 2:    # LH2 fuel cell
    n_engh = 0.6*0.97*0.995*0.85*0.95       # Efficiency engine - hydrogen fuel cell w/o kerosene
    n_engk = 0

#Ratio of hydrogen-to-kerosene
ratio_h = m_fh/(m_fh+m_fk)
ratio_k = m_fk/(m_fh+m_fk)

#PSFC = 0.48*(0.45/(745*3600))
# n_engk = (1/e_fk)*(1/PSFC)
# n_engh = (1/e_fh)*(1/PSFC)

#Aerodynamic charecteristics
LD_crs = 16.2                            # Lift to drag ratio during cruise
L_D = 16                                  #Lift to drag ratio

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 185200   # in m ---> 100nmi
t_E = 45 * 60    # in seconds - endurance time


# #Point A
ranges2 = [0]
plmasses2 = [m_plmax]

#Point B (Range at Max Payload)

m_fB = m_mto - m_oe - m_plmax
m_fhB = ratio_h *m_fB
m_fkB  = ratio_k*m_fB
#R_b1 = n_engh * n_ph * (L_D) * (e_fh/g) * np.log((m_oe + m_plmax + m_fB)/(m_oe + m_plmax))
R_b2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_fhB)/(m_oe + m_pldes))) +  (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe + m_pldes + m_fkB)/(m_oe + m_pldes))))

R_lostB2 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eqB2 = ((R_b2 + R_lostB2)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
R_auxB2 = R_eqB2 - R_b2
R_b = R_b2 - R_auxB2

ranges2 = np.append(ranges2, [R_b])
plmasses2 = np.append(plmasses2, [m_plmax])

#print(R_lostB)
#print(R_eqB)
#print(R_auxB)
#print(m_f)

#Point C (Design Range)
R_lost2 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eq2 = ((R_nom + R_lost2)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
R_aux2 = R_eq2 - R_nom
# m_fh = ((m_mto * (1 - np.exp(-R_nom / (n_engh * n_ph * (e_fh /g) * (L_D)))))*ratio_h)
# m_fk =  ((m_mto * (1 - np.exp(-R_nom / (n_engk * n_pk * (e_fk /g) * (L_D)))))*ratio_k)
R_c2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe + m_pldes + m_fh)/(m_oe + m_pldes))) +  (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe + m_pldes + m_fk)/(m_oe + m_pldes))))   - R_aux2

ranges2 = np.append(ranges2, [R_c2])
plmasses2 = np.append(plmasses2, [m_pldes])

# print(R_lost2)
# print(R_eq2)
# print(R_aux2)
# print("Rc", R_c2)

#Point D (Ferry Range)
R_d2 = ((n_engh * n_ph * (L_D) * (e_fh /g) * np.log((m_oe  + m_fh)/(m_oe )))+ (n_engk * n_pk * (L_D) * (e_fk /g) * np.log((m_oe  + m_fk)/(m_oe )))) - R_aux2

ranges2 = np.append(ranges2, [R_d2])
plmasses2 = np.append(plmasses2, [0])

#print(R_d2)

#Constructing the actual plot

# plotting the points
plt.plot(ranges2, plmasses2, color='red', linewidth=2,
         marker='o', markerfacecolor='red', markersize=5)

# naming the x axis
plt.xlabel('Range (m)')
# naming the y axis
plt.ylabel('Payload (Kg)')

# giving a title to my graph
plt.title('Payload-Range diagram for LH2-Kerosene hybrid')

# function to show the plot
plt.show()
