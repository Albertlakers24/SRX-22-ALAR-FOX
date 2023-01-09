import numpy as np
import matplotlib.pyplot as plt

"""
Constants
lh 
"""


Ch_delta = 1
Cm_delta = 1
Cm_delta_e = -1         #always negative!
xnfree = 1
MACe =1
Ch_delta_t = 1
delta_te =1
delta_te0 =1
fraction = 1
Se = 1
V_cruise= 275*0.51444
V_tailload = np.arange(0,275*0.514444)
V_controlforce = np.arange(0,275*0.51444)
alpha_trim = np.arange(0,10, 0.01)

##Dummy Variables
rho = 1.224             #kg/m^3
MTOW =19200*9.81        #N
xcg = 10.8              #m or 11.3m
xw =0.4*25.6            #m 0.4*l_f

##Dummy Variable for the trim curve
V_trim = np.arange(0,275*0.51444)
Cmac = -0.12            #
alpha_0 = 1*np.pi/180   #radians from Megha
i_h = -2*np.pi/180      #radians to minimize parasite drag
S = 53                  #m^2
MAC = 2.2               #m
lh=11                   #m
Sh=9.1                  #m^2
Vh_V = 1                #for T-tail
Aw =1                   #- ????
Sweep_halfc_w = 22*np.pi/180      #radians ????

##Dummy Variables for the trim curve
M = 0.78                #- Mach number
Ah = 4.584              #-
eta = 0.95              # can be assumed
Sweep_halfc_h = 24*np.pi/180      #radians ????
Cmalpha = 1             #????

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
    CNh_delta = -0.04     #assumed from the graph (slide 19)
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

####Continue hereeeeeee!!!!!!!
print("Slope elevator deflection vs angle =",(delta_e_alpha(alpha_trim[5]-delta_e_alpha(alpha_trim[4])))/(alpha_trim[5]-alpha_trim[4]))



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





#plt.plot(V_tailload, TailLoad(V=V_tailload))
#plt.grid(True)
#plt.xlabel("Velocity (m/s)")
#plt.ylabel("Tail Load (N)")
#plt.title("Tail Loading Graph")
#plt.show()

plt.plot(V_trim, delta_e(V=V_trim))
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Elevator Deflection (m)")
plt.title("Elevator Deflection Curve - Stick Fixed")
plt.show()

plt.plot(alpha_trim, delta_e_alpha(alpha=alpha_trim*np.pi/180))
plt.grid(True)
plt.xlabel("Angle of Attack (m/s)")
plt.ylabel("Elevator Deflection (m)")
plt.title("Elevator Deflection Curve - Stick Fixed")
plt.show()

plt.plot(V_controlforce, ControlForce(V=V_controlforce))
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Control Force (N)")
plt.title("Elevator Control Force Curve")
#plt.show()

plt.plot(V_controlforce, ControlForce(V=V_controlforce))
plt.grid(True)
plt.xlabel("Velocity (m/s)")
plt.ylabel("Control Force (N)")
plt.title("Elevator Control Force Curve")
#plt.show()