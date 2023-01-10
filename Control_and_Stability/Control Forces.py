import numpy as np
import matplotlib.pyplot as plt
from Constants import *
from Class_I_Weight_Estimation.Class_I_weight_estimation_Fuelcell_FINAL import m_mto
from Initial_Aircraft_Sizing.Empennage_Design import l_h, Sh, Ah
from Initial_Aircraft_Sizing.Wing_planform import Sw, c_mac, M_cruise
from Initial_Aircraft_Sizing.cg import x_wing, x_cg_LEMAC

print("-------------------RESULTS FROM CONTROL/STABILITY------------------")

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
M = M_cruise                                            #m/s
Sweep_halfc_h = 24*np.pi/180                            #radians -> GABRIEL
Sweep_halfc_w = 22*np.pi/180                            #radians -> GABRIEL
alpha_0 = 2*np.pi/180                                   #radians -> MEGHA

##Variables
Vh_V = 1                                                #- (for T-tail)
eta = 0.95                                              #-
CNh_delta = 0.04                                        #deg^-1 assumed from the graph (slide 19)

##Dummy Variables
Cmac = 0.12             #-
i_h = 2*np.pi/180       #radians to minimize parasite drag
Cmalpha =-0.05          #deg^-1
xnfix = 1               #m

##Graph Arange velocities and angle of attack
V_tailload = np.arange(V_cruise-V_cruise*0.7,V_cruise+V_cruise*0.4)         #m/s
V_controlforce = np.arange(V_cruise-V_cruise*0.7,V_cruise+V_cruise*0.4)     #m/s
V_trim = np.arange(V_cruise-V_cruise*0.7,V_cruise+V_cruise*0.4)             #m/s
alpha_trim = np.arange(0,10, 0.01)                                          #degrees


##Dummy Variables ControlForce calculations
fraction = 2.25     #rad m^-1
Se = 3              #m^2
MACe = 0.4          #m
Ch_delta = -0.4     #rad^-1
xnfree = 0.483*MAC  #m
delta_te = 0.18127  #rad            have to be re-iterated!
delta_te0 =0        #rad
Ch_delta_t = -0.125 #rad^-1
Vtrim = 100         #m/s

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
    :param Cmac: (-)
    :param CN_alpha: for the horizontal tail plane (rad^-1)
    :param alpha_0: when CL=0 (rad)
    :param i_h: incidence angle to minimize the parasite drag (rad)
    :param Vh_V: fixed at 1 for T-tail (-)
    :return: Cm0 (always positive) (-)
    """
    Cm0 = Cmac - CN_alpha(A = Ah, Sweep_halfc=Sweep_halfc_h)*(alpha_0-i_h)*(Vh_V**2)*(Sh*lh/(S*MAC))
    return Cm0

print("Cm0=", Cm0())
if Cm0() <0:
    print("Cm0 is negative!! Should be positive!")
print("CN_alpha_h=", CN_alpha(A = Ah, Sweep_halfc=Sweep_halfc_h), "used for Cm0")
print("CN_alpha_w", CN_alpha(A = Aw, Sweep_halfc=Sweep_halfc_w), "used for delta_e")

def Cmdelta_e():
    """
    :param CNh_delta: (deg^-1)
    :param Sh and S: surface area (m^2)
    :param MAC: mean aerodynamic chord (m)
    :return: Cmdelta_e: (deg^-1)
    """
    Cmdelta_e = -CNh_delta*(Vh_V**2)*(Sh*lh/(S*MAC))
    return Cmdelta_e

print("Cmdelta_e=", Cmdelta_e())
if Cmdelta_e() >0:
    print("Cmdelta_e is positive! Should be negative!!")

def delta_e(V):
    """
    :param Cm_delta_e:(deg^-1)
    :param Cm0: (-)
    :param Cm_alpha: (rad^-1)
    :param C_N_alpha (rad^-1)
    :param MTOW: maximum take off weight (N)
    :param rho: density (kg/m^3)
    :param S: wing surface area (m^2)
    :param V: velocity (m/s)

    :return: elevator deflection (degrees)
    """
    delta_e = -(1/Cmdelta_e())*(Cm0() + (Cmalpha/CN_alpha(A = Aw, Sweep_halfc=Sweep_halfc_w))*(MTOW/(0.5*S*rho*V**2)))
    return delta_e

print("fraction:", MTOW/(0.5*S*rho*V_cruise**2))

def delta_e_alpha(alpha):
    """
    :param Cmdelta_e: (deg^-1)
    :param Cm0: (-)
    :param Cmalpha: (deg^-1)
    :param alpha: angle of attack (deg)
    :return:
    """
    alpha_0_deg = alpha_0*(1/(np.pi/180))
    delta_e = -(1/Cmdelta_e())*(Cm0() + Cmalpha*(alpha-alpha_0_deg))
    return delta_e

def TailLoad(V,xcg):
    """
    CHECKED
    Tail Load Diagram, used for rotational and vertical stability,
    important for the structural department
    :param Cmac: (m)
    :param rho: density (kg/m^3)
    :param V: velocity (m/s)
    :param MAC: mean aerodynamic chord (m)
    :param MTOW: maximum take off weight (N)
    :param xcg and xw: locations (m)
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

    delta_te0 = ()/Ch_delta_t

    F_velocity_independent = (MTOW/S)*(Ch_delta/Cmdelta_e())*((xcg-xnfree)/MAC)
    F_velocity_dependent= 0.5*rho*V**2*Ch_delta_t*(delta_te-delta_te0)
    a = fraction*Se*MACe*(Vh_V)**2

    Fe = a*(F_velocity_independent - F_velocity_dependent)
    return Fe

def deriv_controlForce():
    deriv_trim = -2*fraction*Se*MACe*(Vh_V**2)*(MTOW/S)*(Ch_delta/Cmdelta_e())*((xcg-xnfree)/MAC)*(1/Vtrim)
    return deriv_trim

print("------------IMPORTANT OUTPUTS FOR STABILITY-----------")
print("Stick Fixed Elevator Deflection")

print("Slope elevator deflection vs angle =",(delta_e_alpha(alpha_trim[5])-delta_e_alpha(alpha_trim[4]))/(alpha_trim[5]-alpha_trim[4]), "deg/deg")
if (delta_e_alpha(alpha_trim[5])-delta_e_alpha(alpha_trim[4]))/(alpha_trim[5]-alpha_trim[4]) >0:
    print("Unstability in the deflection of the elevator due to change in angle of attack :(")
else:
    print("Stability of elevator deflection due to angle of attack :)")

print("Check if close:", "Vtrim=",V_trim[100], "V_cruise=", V_cruise)
print("Slope elevator deflection vs velocity =",(delta_e(V_trim[98])-delta_e(V_trim[97]))/(V_trim[98]-V_trim[97]), "deg/m/s")
if (delta_e(V_trim[98])-delta_e(V_trim[97]))/(V_trim[98]-V_trim[97])<0:
    print("Unstability in the deflection of the elevator due to change in velocity :(")
else:
    print("Stability of elevator deflection due to velocity :)")

print("Check Control Force Stability:")
if deriv_controlForce() <0:
    print("Instability in control force :(")
else:
    print("Stability for Control force derivative at trim speed = ", deriv_controlForce(), "N/m/s")


plt.plot(V_trim, delta_e(V=V_trim))
plt.axvline(x=V_cruise, color="black")
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Elevator Deflection (radians)")
plt.title("Elevator Deflection Curve - Stick Fixed")
plt.legend(["Elevator Deflection", "Vcruise"])
plt.show()

plt.plot(alpha_trim, delta_e_alpha(alpha=alpha_trim*np.pi/180))
plt.grid(True)
plt.xlabel("Angle of Attack (degrees)")
plt.ylabel("Elevator Deflection (radians)")
plt.title("Elevator Deflection Curve - Stick Fixed")
plt.show()

plt.plot(V_tailload, TailLoad(V = V_tailload, xcg=xcg+0.6))
plt.plot(V_tailload, TailLoad(V = V_tailload, xcg=xcg))
plt.axvline(x=V_cruise, color="black")
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Tail Load (N)")
plt.title("Tail Load curve for different velocities and cg locations")
plt.legend(["xcg>xw", "xcg<xw", "Vcruise"])
plt.show()

plt.plot(V_controlforce, ControlForce(V=V_controlforce))
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Control Force (N)")
plt.title("Elevator Control Force Curve - Stick Free")
plt.show()

