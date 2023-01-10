import numpy as np
import pyplot.matplotlib as plt


# Stability Code
lh
c_bar
Vh_V
#Calculation of C_L_alpha_h
A_h
beta
eta
lambda
lambda_0.25 = numpy.radians()
lambda_0.5_c = np.arctan(np.tan(lambda_0.25 - ((4/A)*(0.25*((1-taperratio)/(1 + taperratio))))))
lambda_0.5_c_h = lambda_0.5_c + np.radians(5)
taperratio


C_L_alpha_h = 2* np.pi * A_h / (2 + np.sqrt(4 + ((A_h * beta / eta)*(1+ ((np.tan(lambda_0.5_c_h))**2/beta**2)))))


#Calculation of C_L_alpha_(A-h)

S
Snet
bf
C_L_alpha_w
b


C_L_alpha_(A-h) = (C_L_alpha_w * (1 + (2.15 * (bf/b)))*(Snet/S)) +((np.pi/2)*(bf**2 / S))


#Calculation of Wing Downwash gradient
r
m_tv
K_EA = ((0.1124 + 0.1265*lambda + 0.1766*lambda**2)/r**2) + 0.1024/r +2
K_EA0 = (0.1124/r**2) + (0.1024/r) +2

#Simplification of parts for downwash equation
a = (r**2 + 0.6319 + m_tv**2)**0.5
c = 0.4876*r / (r**2 + m_tv**2)*a

cc = r**2 + 0.7915 + 5.0734*m_tv*82
dd = (1 + (r**2 / cc)**0.3113)*(1 - (1 + mtv**2)**0.5)

downwashgrad_w = ((K_EA * C_L_alpha_w)/(K_EA0 * np.pi * A))*(c + dd)


#Calculation of the aerodynamic center for the aircraft-less configuration

x_ac_w
h_f
b_f
l_f
c_g

x_ac_wf = x_ac_w - ((1.8*b_f * h_f * l_f)/(C_L_alpha_(A-h) * S * c_bar)) + ((0.273*b_f* c_g * (b- b_f))/((1+ taperratio)*c_bar **2 * (b + 2.15*b_f)))*np.tan(lambda_0.25 * np.pi /180)

x_ac_n = kn * 4*(bn * ln)/ S * c_bar * C_L_alpha_(A-h)    #This assumes all 4 propellers are in the same distance from the nose tip
x_ac = x_ac_wf + x_ac_n

#Here goes the final equation for stability, presented in the form y = m_s x + c_s

m_s = 1/((C_L_alpha_h /C_L_alpha_(A-h))* (1-downwashgrad_w)* l_h/c_bar * 1)

c_s = (x_ac -0.05)*m_s

# C O N T R O L L A B I L I T Y


#Calculation of C_L_h

C_L_h = -0.8  #this assumes we will use an adjustable tail (As most commercial airliners do)

#Calculation of C_L_(A-h)
# This can be approached using the CL for landing conditions
W_landing
rho
V_landing

C_L_(A-h) = 2*W_landing / (rho * V_landing **2 * S)


#Calculation of C_m_ac

nac_cont =      #Aproximate from similar aircraft (previous literature/studies)

flap_cont =      #There is an equation but may also be obtained from other aircraft

fus_cont = -1.8 * (1 - 2.5*b_f/l_f)*(np.pi * b_f * h_f * l_f * C_L_0)/(4*S*c_bar * C_L_alpha_(A-h))

C_m_acw = Cm_0airfoil * ((A*np.cos(lambda)**2)/(A + 2*np.cos(lambda)))


C_m_ac = C_m_acw + flap_cont + fus_cont + nac_cont

#Here goes the final equation for Controllability, presented in the form y = m_c x + c_c

m_s = 1/((C_L_h /C_L_(A-h))* l_h/c_bar * 1)

c_s = ((c_m_ac / C_L_(A-h)) - x_ac)*m_s



#Plotting the curves




import matplotlib.pyplot as plt
import numpy as np
x = np.linspace(-1,1,100)   #Vary this as you see
ys = m_s * x + c_s
plt.plot(x, ys, '-r', label='Stability')
yc = m_c * x + c_c
plt.plot(x, yc, 'g', label = 'Controllability')
plt.title('Scissor Plot')
plt.xlabel('x', color='#1C2833')
plt.ylabel('y', color='#1C2833')
plt.legend(loc='upper left')
plt.grid()
plt.show()