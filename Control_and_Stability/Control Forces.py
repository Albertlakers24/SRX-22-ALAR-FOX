import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Initial_Aircraft_Sizing.Empennage_Design import l_h, Sh, Ah
from Initial_Aircraft_Sizing.Wing_planform import Sw, c_mac
from Initial_Aircraft_Sizing.cg import x_wing, x_cg_LEMAC
"""
Imported Variables come from respective files
"""



##Imported Variables
lh= l_h                                                 #m
rho = ISA_calculator(h = h_cruise, dt=dt_cruise)[2]     #kg/m^3
S = Sw                                                  #m^2
MAC = c_mac                                             #m
MTOW =m_mto*g                                           #N
xcg = x_cg_LEMAC                                        #m
xw =x_wing                                              #m
Sh=Sh                                                   #m^2
Ah = Ah                                                 #-
Aw = A                                                  #-
M = V_cruise/a                                          #m/s
Sweep_halfc_h = 24*np.pi/180                            #radians -> GABRIEL
Sweep_halfc_w = 22*np.pi/180                            #radians -> GABRIEL
alpha_0 = 1*np.pi/180                                   #radians -> MEGHA

##Dummy Variables
Cmac = -0.12            #-
CNh_delta = -0.04       #- assumed from the graph (slide 19)
Vh_V = 1                #for T-tail
i_h = -2*np.pi/180      #radians to minimize parasite drag
eta = 0.95              # can be assumed
Cmalpha =1              #
xnfix = 1               #




##Graph Arange velocities and angle of attack
V_tailload = np.arange(0,275*0.514444)
V_controlforce = np.arange(0,275*0.51444)
V_trim = np.arange(0,275*0.51444)
alpha_trim = np.arange(0,10, 0.01)





def CN_alpha(A,Sweep_halfc):
    """
    (CHECKED)
    Equation for the lift gradients, as CN_alpha can be assumed equal to CL_alpha
    :param A: Aspect Ratio (-)
    :param Sweep_halfc: sweep at half chord (radians)
    :param M: Mach number (-)
    :param eta: airfoil efficiency coefficient (always 0.95)
    :return: C_N_alpha (rad^-1)
    """
    beta = np.sqrt(1-M**2)
    CN_alpha = (2*np.pi*A)/(2+np.sqrt(4+((A*beta/eta)**2)*(1+((np.tan(Sweep_halfc))**2)/beta**2)))
    return CN_alpha

def Cm0():
    """
    :param Cmac:
    :param CN_alpha: for the horizontal tail plane
    :param alpha_0: when CL=0
    :param i_h: incidence angle to minimize the parasite drag
    :param Vh_V: fixed at 1 for T-tail
    :return: Cm0 (always positive)
    """
    Cm0 = Cmac - CN_alpha(A = Ah, Sweep_halfc=Sweep_halfc_h)*(alpha_0+i_h)*(Vh_V**2)*(Sh*lh/S*MAC)
    return Cm0

print("Cm0=", Cm0())
print("CN_h=", CN_alpha(A = Ah, Sweep_halfc=Sweep_halfc_h))

if Cm0() <0:
    print("Cm0 is negative!! Should be positive!")

def Cmdelta_e():
    Cmdelta_e = -CNh_delta*(Vh_V**2)*(Sh*lh/(S*MAC))
    return Cmdelta_e

print("Cmdelta_e=", Cmdelta_e())

def delta_e(V):
    """
    :param V: velocity (m/s)
    :param Cm_delta_e:XXXXXX
    :param Cm0: XXXXXXX
    :param Cm_alpha: XXXXXX
    :param MTOW: maximum take off weight (N)
    :param rho: density (kg/m^3)
    :param S: wing surface area (m^2)

    :return: elevator deflection (m)
    """
    delta_e = -(1/Cmdelta_e())*(Cm0() + (Cmalpha/CN_alpha(A = Aw, Sweep_halfc=Sweep_halfc_w))*(MTOW/(0.5*rho*V**2*S)))
    return delta_e

def delta_e_alpha(alpha):
    delta_e = -(1/Cmdelta_e())*(Cm0() + Cmalpha*(alpha-alpha_0))
    return delta_e


print("------------IMPORTANT OUTPUTS-----------")
print("Stick Fixed Elevator Deflection")

print("Slope elevator deflection vs angle =",(delta_e_alpha(alpha_trim[5])-delta_e_alpha(alpha_trim[4]))/(alpha_trim[5]-alpha_trim[4]))



def TailLoad(V):
    """
    Tail Load Diagram, used for rotational and vertical stability,
    important for the structural department
    :param V: Range from 0 to XXXX (m/s)
    :return: Tail Load (N)
    """
    Nh = (1/lh)*(Cmac*0.5*rho*V**2*S*MAC + MTOW*(xcg-xw))
    return Nh

def CNh():
    CNh = TailLoad(V = V_cruise)/(0.5*rho*V_cruise**2*S)
    return CNh

def ControlForce(V):
    """
    Control Curve, used for elevator control force stability (dFe/dV)Fe=0 >0
    :param V: range from 0 to XXXX (m/s)
    :return: Tail Load (N)
    """
    F_velocity_independent = (MTOW/S)*(Ch_delta/Cm_delta_e)*((xcg-xnfree)/MAC)
    F_velocity_dependent= 0.5*rho*V**2*Ch_delta_t*(delta_te-delta_te0)
    a = fraction*Se*MACe*(Vh_V)**2

    Fe = a*(F_velocity_independent - F_velocity_dependent)
    return Fe


plt.plot(V_trim, delta_e(V=V_trim))
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Elevator Deflection (m)")
plt.title("Elevator Deflection Curve - Stick Fixed")
#plt.show()

plt.plot(alpha_trim, delta_e_alpha(alpha=alpha_trim*np.pi/180))
plt.grid(True)
plt.xlabel("Angle of Attack (m/s)")
plt.ylabel("Elevator Deflection (m)")
plt.title("Elevator Deflection Curve - Stick Fixed")
#plt.show()

#plt.plot(V_controlforce, ControlForce(V=V_controlforce))
#plt.grid(True)
#plt.xlabel("Velocity (m/s)")
#plt.ylabel("Control Force (N)")
#plt.title("Elevator Control Force Curve")
#plt.show()

