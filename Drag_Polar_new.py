import numpy as np
import matplotlib.pyplot as plt
# from Class_I_Weight_Estimation import Cd0, e, A
#prop_type = 1   # 1 for LH2 (both combustion & fuel cell), 2 for Electric Hybrid (Both configurations)

# INPUT CL_Max for each flight phase  ---> Call from Class I ---- EVENTUALLY
CL_max_take0ff = 2.1
CL_max_cruise = 1.9
CL_max_landing = 2.6

# INPUTS from Class I
##Cdo calculations
Psi = 0.0075                    #Parasite drag dependent on the lift coefficient (value based on Roelof reader p.46)
phi = 0.97                      #span efficiency factor (value based on Roelof reader p.46)
A = 12                           #Aspect Ratio (12-14) #Reference to ATR 72
e_cruise = 1/(np.pi*A*Psi+(1/phi))
Cfe = 0.0030                     #equivalent skin friction coefficient -> depending on aircraft from empirical estimation
Swet_S = 6.1                     #(6.0-6.2) wetted area ratios -> depending on airframe structure
CD0_Cruise = Cfe * Swet_S

# FLAP DEFLECTION --> FROM ADSEE READER
delta_f_take_off = 15       # TAKE Flap deflection - in degrees
delta_f_landing = 35        # LANDING Flap deflection - in degrees

# Make CL Lists for each flight phase
step = 0.01
CL_list_takeOff = np.arange(0, CL_max_take0ff + step, step)
CL_list_cruise = np.arange(0, CL_max_cruise + step, step)
CL_list_landing = np.arange(0, CL_max_landing + step, step)

# Make CD Lists for each flight phase
CD_list_cruise = []
#LH2
CD_list_takeOff_woLG_LH2 = []
CD_list_takeOf_withLG_LH2 = []
CD_list_landing_LH2 = []
# ELECTRIC HYBRID
CD_list_takeOff_woLG_EH = []
CD_list_takeOf_withLG_EH = []
CD_list_landing_EH = []

# CDO Change due to flap deflection
CD0_Land_flapDef = (13 *10**(-4) * delta_f_landing)
CD0_TO_flapDef = (13 *10**(-4) * delta_f_take_off)

# CD0 change due to Landing gear
CD0_LG = 175*10**(-4)

# CD0 change per configuration
Change_CD0_TO_woLG = CD0_TO_flapDef
Change_CD0_TO_WithLG = CD0_TO_flapDef + CD0_LG
Change_CD0_landing = CD0_Land_flapDef + CD0_LG

# CDO per configuration
CD0_TO_woLG = CD0_Cruise + Change_CD0_TO_woLG
CD0_TO_WithLG = CD0_Cruise + Change_CD0_TO_WithLG
CD0_Landing = CD0_Cruise + Change_CD0_landing

# Updating CDO list
# CD_list_takeOff_woLG.append(CD0_TO_woLG)
# CD_list_takeOf_withLG.append(CD0_TO_WithLG)
# CD_list_cruise.append(CD0_Cruise)
# CD_list_landing.append(CD0_Landing)

# Change in 'e' per configuration -- for flap deflection
#if prop_type == 1:      # LH2 propulsion - ONLY WING MOUNTED
Change_e_LH2_TO = 0.0026 * delta_f_take_off
Change_e_LH2_Land = 0.0026 * delta_f_landing
e_LH2_TO = e_cruise + Change_e_LH2_TO
e_LH2_Land = e_cruise + Change_e_LH2_Land

#if prop_type == 2:      # Electric Hybrid propulsion - WING + FUSELAGE MOUNTED
Change_e_EH_TO = (0.0026 * delta_f_take_off) + (0.0046 * delta_f_take_off)
Change_e_EH_Land = (0.0026 * delta_f_landing) + (0.0046 * delta_f_take_off)
e_EH_TO = e_cruise + Change_e_EH_TO
e_EH_Land = e_cruise + Change_e_EH_Land

# DRAG LH2 PROPULSION - ONLY WING MOUNTED!!

red_CD = 0.9
for i in range(len(CL_list_cruise)):
    CDi_cruise = (CL_list_cruise[i])**2/(np.pi *A*e_cruise)
    CD_cruise = (CD0_Cruise + CDi_cruise) * red_CD
    CD_list_cruise.append(CD_cruise)
CL_CD = CL_list_cruise / CD_list_cruise
max_CL_CD = max(CL_CD)
print(CL_CD)
print(max_CL_CD, "max L/D")
index = np.where(CL_CD == max_CL_CD)
print(CL_list_cruise[index])
# print(CL_CD.index(max_CL_CD))


# LH2 CASES

for i in range(len(CL_list_takeOff)):
    CDi_TO_LH2 = (CL_list_takeOff[i])**2/(np.pi *A*e_LH2_TO)
    CD_TO_woLG = CD0_TO_woLG + CDi_TO_LH2
    CD_list_takeOff_woLG_LH2.append(CD_TO_woLG)

for i in range(len(CL_list_takeOff)):
    CDi_TO_LH2 = (CL_list_takeOff[i])**2/(np.pi *A*e_LH2_TO)
    CD_TO_WithLG = CD0_TO_WithLG + CDi_TO_LH2
    CD_list_takeOf_withLG_LH2.append(CD_TO_WithLG)

for i in range(len(CL_list_landing)):
    CDi_land_LH2 = (CL_list_landing[i])**2/(np.pi *A*e_LH2_Land)
    CD_Land_LH2 = CD0_Landing + CDi_land_LH2
    CD_list_landing_LH2.append(CD_Land_LH2)

# ELECTRIC HYBRID CASES

for i in range(len(CL_list_takeOff)):
    CDi_TO_EH = (CL_list_takeOff[i]) ** 2 / (np.pi * A * e_EH_TO)
    CD_EH_woLG = CD0_TO_woLG + CDi_TO_EH
    CD_list_takeOff_woLG_EH.append(CD_EH_woLG)

for i in range(len(CL_list_takeOff)):
    CDi_TO_EH = (CL_list_takeOff[i]) ** 2 / (np.pi * A * e_EH_TO)
    CD_TO_WithLG = CD0_TO_WithLG + CDi_TO_EH
    CD_list_takeOf_withLG_EH.append(CD_TO_WithLG)

for i in range(len(CL_list_landing)):
    CDi_land_EH = (CL_list_landing[i]) ** 2 / (np.pi * A * e_EH_Land)
    CD_Land_EH = CD0_Landing + CDi_land_EH
    CD_list_landing_EH.append(CD_Land_EH)

'''plt.plot(CD_list_cruise, CL_list_cruise, color = 'red', label = 'Cruise')
plt.plot(CD_list_takeOff_woLG_EH, CL_list_takeOff, color = 'green', label = 'Take-off, Landing gear up')
plt.plot(CD_list_takeOf_withLG_EH, CL_list_takeOff, color = 'orange', label = 'Take-off, Landing gear down')
plt.plot(CD_list_landing_EH, CL_list_landing, color = 'blue', label = 'Landing')
plt.legend()
plt.xlabel('Drag co-efficient $C_{D}$ [-]')   # naming the x axis
plt.ylabel('Lift co-efficient $C_{L}$ [-]')  # naming the y axis
plt.title('Drag polar for Electric Hybrid propulsion') # giving a title to my graph
plt.show()'''

plt.plot(CD_list_cruise, CL_list_cruise, color = 'red', label = 'Cruise')
plt.plot(CD_list_takeOff_woLG_LH2, CL_list_takeOff, color = 'green', label = 'Take-off, Landing gear up')
plt.plot(CD_list_takeOf_withLG_LH2, CL_list_takeOff, color = 'orange', label = 'Take-off, Landing gear down')
plt.plot(CD_list_landing_LH2, CL_list_landing, color = 'blue', label = 'Landing')
plt.legend()
plt.xlabel('Drag co-efficient $C_{D}$ [-]')   # naming the x axis
plt.ylabel('Lift co-efficient $C_{L}$ [-]')  # naming the y axis
plt.title('Drag polar for hydrogen combustion') # giving a title to my graph
plt.show()
