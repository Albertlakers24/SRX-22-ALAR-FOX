import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Wing_Loading_Diagram import V_approach_stall, beta_s_land_fc
from Class_I_Weight_Estimation.Class_I_fuel_cell import m_MTOW
from Initial_Aircraft_Sizing.Empennage_Design import l_h, Ah
from Initial_Aircraft_Sizing.Wing_planform import c_mac, M_cruise, Sw, A, b, c_r, taper
from Control_and_Stability.Control_Forces import Sweep_halfc_h, Sweep_halfc_w
from Initial_Aircraft_Sizing.Fuselage import D_outer, l_f

##Imported Variable
lh = l_h
c_bar = c_mac
Vh_V = 1
M = M_cruise   #m/s
A_h = Ah
A_w = A
beta = np.sqrt(1-M**2)
eta = 0.95
lambdahalf_h = 2.05*np.pi/180                          #radians -> GABRIEL
lambdahalf_w = 1.81*np.pi/180                          #radians -> GABRIEL
S = Sw
bf = D_outer
b = b
cr = c_r
hf = D_outer             #m
cg = Sw/b                #m mean geometric chord
lambda_quarterchord =0   #rad
SweepLE = 1.81                #degrees at the LE -> GABRIEL
rho = ISA_calculator(h = h_cruise, dt=dt_cruise)[2]     #kg/m^3
W_landing = m_MTOW*beta_s_land_fc
C_L_h_adj = -0.8             #this assumes we will use an adjustable tail (As most commercial airliners do)
C_L_h_mov = -1              #for a full moving tail
C_L_h_fix = -0.35*A_h**(1/3)    #fixed


##Dummy Variables
vt = 5                   #vertical distance between ac of wing and tail
ln = -10                   #?? m distance of front nacelle to quarter chord
bn = 3                   #?? m width of nacelle
l_fn = 10.24             #?? m TBD nose to nacelle ?? - i DONT TRUST GABRIEL1!!!!
V_landing = V_approach      #or use V_approach_stall
print("landing weight",W_landing)
SM = 0.05               #normal value
xcg_MAC_gear = 0.4      #??

def C_L_alpha(A, lambdahalf):
    C_L_alpha = 2 * np.pi * A / (2 + np.sqrt(4 + ((A * beta / eta) * (1 + ((np.tan(lambdahalf)) ** 2 / beta ** 2)))))
    return C_L_alpha

def CL_alpha_Ah():
    Snet = S - cr*bf
    C_L_alpha_Ah = (C_L_alpha(A= A, lambdahalf=lambdahalf_w) * (1 + (2.15 * (bf / b))) * (Snet / S)) + ((np.pi / 2) * (bf ** 2 / S))
    return C_L_alpha_Ah

def Downwash():
    r = lh / (b / 2)
    mtv = 2*(vt/b)
    K_EA = ((0.1124 + 0.1265 * SweepLE + 0.1766 * SweepLE ** 2) / r ** 2) + 0.1024 / r + 2
    K_EA0 = (0.1124 / r ** 2) + (0.1024 / r) + 2
    downwash = (K_EA/K_EA0)*(r/(r**2+mtv**2) * (0.4876/np.sqrt(r**2+0.6319+mtv**2)) +
    (1+ (r**2/(r**2+0.7915+5.073*mtv**2))**0.3113)*(1-np.sqrt(mtv**2/(1+mtv**2))))*(C_L_alpha(A= A_w, lambdahalf=lambdahalf_w)/(np.pi*A_w))
    return downwash

def AC_location():
    x_ac_w = 0.26  # DOROTHEA GRAPH
    x_ac_f1 = -((1.8 * bf * hf * l_fn) / (CL_alpha_Ah() * S * c_bar))
    x_ac_f2 = ((0.273 * bf * cg * (b - bf)) / ((1 + taper) * c_bar ** 2 * (b + 2.15 * bf))) * np.tan(lambda_quarterchord)
    x_ac_n = -4*(bn**2 * ln)/ (S * c_bar * CL_alpha_Ah())
    x_ac = x_ac_w+x_ac_f1+x_ac_f2+x_ac_n
    return x_ac

def C_m_AC():
    """
    To be done: flap, fus and nac contributions!!
    :return:
    """
    Cm_0airfoil = 0.1
    C_L_0 = 0.1
    C_m_acw = Cm_0airfoil* (A*(np.cos(SweepLE)**2))/(A+2*np.cos(SweepLE))
    flap_cont =  1                   #There is an equation but may also be obtained from other aircraft
    fus_cont = -1.8 * (1 - 2.5*bf/l_f)*(np.pi * bf * hf * l_f * C_L_0)/(4*S*c_bar * CL_alpha_Ah())
    nac_cont = 1                      #Aproximate from similar aircraft (previous literature/studies)
    C_m_ac = C_m_acw + flap_cont + fus_cont + nac_cont
    return C_m_ac

def C_L_Ah():
    C_L_Ah = 2*W_landing/(rho*S*V_landing**2)
    return C_L_Ah

x = np.arange(-0.5,0.8,0.01)             ##Note should be xcg/MAC

#Final equation for STABILITY, presented in the form y = m_s x + c_s
m_s = 1/((C_L_alpha(A = Ah, lambdahalf=lambdahalf_h) /CL_alpha_Ah())* (1-Downwash())* lh/c_bar * (Vh_V**2))
c_s = (AC_location()-0.05)*m_s
c_s_SM = AC_location()*m_s
ys = m_s * x + c_s
ys_SM = m_s*x +c_s_SM

#Here goes the final equation for CONTROLLABILITY, presented in the form y = m_c x + c_c
def y_c(C_L_h):
    m_c = 1 / ((C_L_h / C_L_Ah()) * lh / c_bar * (Vh_V ** 2))
    c_c = ((C_m_AC() / C_L_Ah()) - AC_location()) * m_c
    yc = m_c * x + c_c
    return yc






#cg locations


#Plotting the curves
plt.plot(x, ys, '-r', label='Neutral Stability Line')
plt.plot(x,ys_SM, 'g', label="Stability Line")
plt.plot(x, yc, 'b', label = 'Controllability')
plt.axvline(x=xcg_MAC_gear, color="black")
plt.title('Scissor Plot')
plt.xlabel('xcg/MAC', color='#1C2833')
plt.ylabel('Sh/S', color='#1C2833')
plt.legend(loc='upper left')
plt.xlim([-0.4, 0.8])
plt.ylim([0,2])
plt.grid(True)
plt.show()

#x = np.linspace(-1,1,100)   #Vary this as you


"""
Graph for effect of wing shift on cg travel see slide 30 Lecture 8
"""