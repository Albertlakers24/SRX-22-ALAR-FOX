import numpy as np
from matplotlib import pyplot as plt
n_eng =                                                    #Engine efficiency
n_p                                                        #Propulsive efficiency
L/D                                                        #Lift over drag
e_f                                                        #
g                                                        #
m_oe                                                   #
m_plmax                                                   #
m_f                                                   #
m_mto                                                   #
#Point A
ranges = [0]
plmasses = [m_pl]
#Point B (Range at Max Payload)
R_b = n_eng * n_p * (L/D) * (e_f /g) * np.log((m_oe + m_plmax + m_f)/(m_oe + m_pl))

m_f = m_mto * (1 - np.exp(-R_b / (n_eng * n_p * (e_f /g) * (L/D))))

ranges = np.append(ranges, [R_b])
plmasses = np.append(plmasses, [m_pl])

print(ranges)

#Point C (Design Range)

R_lost = (1/0.7) * (L/D)_crs * (h_cr + ((V_cr **2)/(2*g)))
R_eq = ((R_nom + R_lost)*(1+f_con)) + (1.2*R_div) + (t_E * V_cr)
R_aux = R_eq - R_nom
m_pldes

R_c = (n_eng * n_p * (L/D) * (e_f /g) * np.log((m_oe + m_pl + m_f)/(m_oe + m_pldes))) - R_aux

ranges = np.append(ranges, [R_c])
plmasses = np.append(plmasses, [m_pldes])

#Point D (Ferry Range)


