import numpy as np
from matplotlib import pyplot as plt

prop_system = 2  # put 1 for LH2 fuel cell, 2 for LH2 combustion

g = 9.80665

# Masses Aircraft --> Payload - constant for all configurations

m_pldes = 5443  # Payload mass designed for
m_plmax = m_pldes*1

# Masses Aircraft --> Max, empty, fuel - per propulsion system
# fuel = 555, tank = 666
if prop_system == 1:  # LH2 fuel cell
    m_fh = 276 # LH2 mass in kg
    m_tank = 0.5*m_fh
    m_oe = (12727.3) + m_tank # Operating empty mass + tank mass
    m_mto = (19391 - 555 - 666) + m_fh + m_tank
    m_fk = 0   # Kerosene mass in kg
if prop_system == 2:    # LH2 Combustion
    m_mto = 19223  # in kg --> Max take off
    m_oe = 12512  # in kg --> Operating empty
    m_fh = 1268  # in kg --> Fuel mass
    m_fk = 0
if prop_system == 3:   # LH2_Kerosene fuel cell
    m_oe = 14643  # Operating empty mass
    m_mto = 21207
    m_fk = 390 # Kerosene mass in kg
    m_fh = 445 # LH2 mass in kg
# Propulsion charecteristics
n_pk =  0.85                            #Propulsive efficiency twin turboprop
n_ph = 0.85                             #Propulsive efficiency twin turboprop
e_fk = 42.9 * 10**6                     # Specific Energy - kerosene
e_fh = 120 * 10**6                      # specific energy - hydrogen
print(m_mto)
if prop_system == 1:    # LH2 fuel cell
    n_engh = 0.6*0.97*0.995*0.85*0.95       # Efficiency engine - hydrogen fuel cell w/o kerosene
    n_engk = 0
if prop_system == 2:    # LH2 Combustion
    n_engh = 0.3
    n_engk = 0
if prop_system == 3:    # LH2_Kerosene fuel cell
    n_engk = 0.39*0.9*0.99                  # Efficiency engine - kerosene in fuel cell system
    n_engh = 0.6*0.97*(0.995**2)*0.85*0.95       # Efficiency engine - hydrogen fuel cell w kerosene


#Ratio of hydrogen-to-kerosene
ratio_h = m_fh/(m_fh+m_fk)
ratio_k = m_fk/(m_fh+m_fk)

#PSFC = 0.48*(0.45/(745*3600))
# n_engk = (1/e_fk)*(1/PSFC)
# n_engh = (1/e_fh)*(1/PSFC)

#Aerodynamic charecteristics
LD_crs = 16.7                            # Lift to drag ratio during cruise
#L_D = 16                                  #Lift to drag ratio

# Mission Charecteristics
V_cr = 141.471   # in m/s
h_cr = 8534.4    # in m
R_nom = 1852000/2  # in m ---> Design range - 1000nmi
f_con = 5/100
R_div = 0 #185200   # in m ---> 100nmi
t_E =  0 #45 * 60    # in seconds - endurance time


# #Point A


#Point B (Range at Max Payload)

m_fB = m_mto - m_oe - m_plmax
m_fhB = ratio_h *m_fB
m_fkB  = ratio_k*m_fB
#R_b1 = n_engh * n_ph * (L_D) * (e_fh/g) * np.log((m_oe + m_plmax + m_fB)/(m_oe + m_plmax))
R_b2 = ((n_engh * n_ph * (LD_crs) * (e_fh /g) * np.log((m_oe + m_pldes + m_fhB)/(m_oe + m_pldes))) +  (n_engk * n_pk * (LD_crs) * (e_fk /g) * np.log((m_oe + m_pldes + m_fkB)/(m_oe + m_pldes))))

R_lostB2 = (1/0.7) * (LD_crs) * (h_cr + ((V_cr **2)/(2*g)))
R_eqB2 = ((R_b2 + R_lostB2)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
R_auxB2 = R_eqB2 - R_b2
R_b = R_b2 - R_auxB2
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
R_c = ((n_engh * n_ph * (LD_crs) * (e_fh /g) * np.log((m_oe + m_pldes + m_fh)/(m_oe + m_pldes))) +  (n_engk * n_pk * (LD_crs) * (e_fk /g) * np.log((m_oe + m_pldes + m_fk)/(m_oe + m_pldes))))   - R_aux2


# print(R_lost2)
# print(R_eq2)
# print(R_aux2)
# print("Rc", R_c2)

#Point D (Ferry Range)
R_d = ((n_engh * n_ph * (LD_crs) * (e_fh /g) * np.log((m_oe  + m_fh)/(m_oe )))+ (n_engk * n_pk * (LD_crs) * (e_fk /g) * np.log((m_oe  + m_fk)/(m_oe )))) - R_aux2

#print(R_d2)

ranges = np.array([0, R_b/1852, R_c/1852, R_d/1852])
plmasses = np.array([m_plmax, m_plmax, m_pldes,0])

print(ranges, plmasses)
#Constructing the actual plot

# plotting the points


# giving a title to my graph
if prop_system == 1:
    plt.plot(ranges, plmasses, color='red', linewidth=2,
             marker='o', markerfacecolor='red', markersize=5)
    plt.title('Payload-Range diagram for LH2 fuel cell')
if prop_system == 2:
    plt.plot(ranges, plmasses, color='green', linewidth=2,
             marker='o', markerfacecolor='green', markersize=5)
    plt.title('Payload-Range diagram for LH2 Combustion')
if prop_system == 3:
    plt.plot(ranges, plmasses, color='blue', linewidth=2,
             marker='o', markerfacecolor='blue', markersize=5)
    plt.title('Payload-Range diagram for LH2-Kerosene Fuel cell')

n = ['A','B','C','D']
for i, txt in enumerate(n):
    plt.annotate(txt, (ranges[i], plmasses[i]))

# naming the x axis
plt.xlabel('Range (nmi)')
# naming the y axis
plt.ylabel('Payload (kg)')

# function to show the plot
plt.show()
