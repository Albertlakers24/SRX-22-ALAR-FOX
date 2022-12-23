import numpy as np
import matplotlib.pyplot as plt

"""
Verification Constants:
RF_CO2 is in W/m^2
horizon is the time horizon- 30yrs (Dallara) for short-lived species and 100 yrs for long-lived species
Time span - from 2035 to 2065
massflow is in kg/s
t_year is time flown in a year in seconds
number_aircraft is the expected fleet size
"""

RF_CO2 = 3.7
horizon = 100
time = np.arange(2035, 2035+horizon, 1)
time_horizon = np.arange(0,horizon, 1)
massflow = 0.031153898
t_year = 747900*12
number_aircraft = 100
Distance_year = 600000

"""Temperature Model Parameters from Dallara:
S = 3 in K
tau in yrs"""
S = 3
alpha_t = 0.595
tau_t1 = 8.4
tau_t2 = 409.5


def convolve(signal, response):
    """
    Determines the convolution integral for the given signal and response (or
    filter) function.

    :param signal: Signal, such as emissions over years
    :param response: Response function to be applied to signal, such as decay
    or other filter
    :return: Convolution integral of provided signal and response function
    """
    output = np.convolve(signal, response, 'full')
    return output[:signal.size]


def RF_contrails(RF_CO2):
    """
    :param ratio: Fixed ratio (W/m^2/nmi)                           -DONE maybe changed for FC
    :param lenght: dependent on mission  50,000 (nmi/year)          -DONE
    :param forcing_factor: dependent on altitude (-)                -TO BE REWRITTEN DEPENDENT ON ALTITUDE
    :param RF_CO2: doubling coefficient (W/m^2)                     -DONE
    :return: RF_contrail_norm: normalized RF (W/m^2/year)           -DONE
    """
    ratio = 2.21*10**-12
    length = Distance_year*number_aircraft
    Eff = 0.59
    forcing_factor= 0.83
    t_test = np.ones(horizon)

    RF_contrail = ratio*length*forcing_factor
    RF_contrail_norm = t_test*Eff * RF_contrail / RF_CO2

    return RF_contrail_norm

def RF_water(massflow, t_year, RF_CO2):
    """
    Calculation of the normalized RF for water emissions
    :param ratio: fixed RF parameter for water (W/m^2/kg)           -DONE
    :param EI: Emission index for water (kg/kg)                     -DONE
    :param massflow: engine (kg/s)                                  -DONE
    :param RF_CO2: doubling coefficient (W/m^2)                     -DONE
    :return: normalized RF (W/m^2/year)                             -DONE
    """
    EI = 2.6*1.26
    ratio = 7.43 * 10 ** -15
    Eff = 1.14

    massflow_month = massflow*t_year*number_aircraft
    E = EI*massflow_month

    t_test = np.ones(horizon)
    RF_water = ratio*E
    RF_water_norm = t_test*Eff*RF_water/RF_CO2
    return RF_water_norm

def G_factor(time):
    """
    :param time: time array (years)                                   -DONE
    :return: G-factor (/year)                                         -DONE
    """

    G_factor = S*((alpha_t/tau_t1)*np.exp(-time/tau_t1) + ((1-alpha_t)/tau_t2)*np.exp(-time/tau_t2))
    return G_factor

#G_factor = S*((alpha_t/tau_t1)*np.exp(-time/tau_t1) + ((1-alpha_t)/tau_t2)*np.exp(-time/tau_t2))

def Temperature_response(RF, time):
    """
    Estimates of temperature response
    :param RF is the total RF for contrails and water
    :param time is the time in years
    :return: Temperature response (K/year)
    """
    result = convolve(signal=RF,response=G_factor(time=time_horizon))
    return result

"""RF_total is dependent on time (per years)"""
RF_total = RF_contrails(RF_CO2=RF_CO2) + RF_water(massflow=massflow, t_year=t_year, RF_CO2=RF_CO2)
print("RF_total=", RF_total, "per year")

"""Result on Delta T"""
result_Temperatureresponse = Temperature_response(RF=RF_total,time=time)
print("Temperature Change=",result_Temperatureresponse)


plt.plot(time, result_Temperatureresponse)
plt.grid(True)
plt.xlabel("Year")
plt.ylabel("Average Temperature Change")
plt.title("Average Temperature Response for a time horizon of 30 years for an entire fleet")
plt.show()

plt.plot(time, G_factor(time=time_horizon))
plt.grid(True)
plt.xlabel("Year")
plt.ylabel("Scaled temperature impulse response")
plt.title("Scaled Temperature Impulse Response function for designated time")
plt.show()

#RF_test needs to be made an array
print("length RF=", len(RF_total))
print("length time", len(time))
print("length G=", len(G_factor(time=time)))
print("type RF=", type(RF_total))
print("type time=", type(time))
print("type G_factor=", type(G_factor(time=time)))