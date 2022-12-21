import numpy as np


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


def concentration_response(t_years,
                           alpha_i=np.array(
                               [0.067, 0.1135, 0.152, 0.097, 0.041]),
                           tau_i=np.array([np.inf, 318.8, 79.8, 18.8, 1.7])):
    """
    Delta CO2 concentration response due to CO2 emissions.

    :param t_years: Array with time, typically years
    :type t_years: np.ndarray
    :param alpha_i: Decay factor in ppbv/Tg
    :param tau_i: Lifetime in years
    :return: Delta CO2 concentration response function
    """
    tn = t_years - t_years[0]  # + 1
    alpha_in = alpha_i * 1e-3 * 1e-9  # convert from ppbv/Tg to ppmv/kg
    try:
        one = np.ones((5, tn.size)) * -1.
    except AttributeError:
        one = np.ones((5, 1)) * -1.
    one = one * tn
    one = one.transpose()
    one = one / tau_i
    one = np.exp(one)

    one = one * alpha_in

    one = np.sum(one, axis=1)
    return one


def sausen_concentration_carbon_dioxide(t_years, emissions_co2, co2_0=380.,
                                        alpha_i=np.array(
                                            [0.067, 0.1135, 0.152, 0.097,
                                             0.041]),
                                        tau_i=np.array(
                                            [np.inf, 318.8, 79.8, 18.8, 1.7])):
    """
    Estimates CO2 concentration in atmosphere in ppmv.

    :param t_years: Array with time, typically years
    :type t_years: np.ndarray
    :param emissions_co2: CO2 emissions per year in t_year, in kg
    :type emissions_co2: np.ndarray
    :param co2_0: Pre-industrial or background CO2 concentration in
        ppmv
    :type co2_0: float
    :param alpha_i: Decay factor in ppbv/Tg
    :param tau_i: Lifetime in years
    :return: CO2 concentration in atmosphere in ppmv
    """
    delta = convolve(signal=emissions_co2,
                     response=concentration_response(t_years=t_years,
                                                     alpha_i=alpha_i,
                                                     tau_i=tau_i))
    print("delta=", np.sum(delta), "for many years")

    # delta = np.cumsum(delta)
    return delta + co2_0


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # Data CO2 from Sausen and Schumann (2000)
    t_test = np.arange(1940, 2031, 1)
    e_co2_test = []
    fal_1940_1995 = np.array(
        [7.7, 8.3, 9.0, 9.7, 10.5, 11.3, 12.3, 13.2, 14.3, 15.4,
         16.7, 18.0, 19.4, 21.0, 22.7, 24.5, 26.5, 28.6, 30.9, 33.3,
         36.0, 39.5, 43.5, 45.9, 48.0, 51.3, 55.6, 65.6, 74.3, 77.8,
         78.0, 90.0, 96.0, 99.4, 96.0, 96.1, 96.4, 102.1, 105.7, 110.1,
         110.9, 109.3, 110.5, 112.0, 119.5, 123.4, 129.9, 135.6, 141.4, 146.5,
         146.9, 143.4, 142.0, 144.1, 150.0, 154.3])
    for i, ti in enumerate(t_test):
        if ti <= 1995:
            e_co2_test.append(fal_1940_1995[i])
        elif 1995 < ti <= 2015:
            emi = ((278.7 - 154.3) / (2015 - 1995)) * (ti - 1995) + 154.3
            e_co2_test.append(emi)
        elif 2015 < ti <= 2050:
            emi = ((405.1 - 278.7) / (2050 - 2015)) * (ti - 2015) + 278.7
            e_co2_test.append(emi)
        elif 2050 < ti <= 2100:
            emi = ((666.2 - 405.1) / (2100 - 2050)) * (ti - 2050) + 405.1
            e_co2_test.append(emi)

    e_co2_test = np.array(e_co2_test) * 1e9  # convert to kilograms

    test_co2_concentration = sausen_concentration_carbon_dioxide(
        t_years=t_test,
        emissions_co2=e_co2_test)

    plt.plot(t_test, test_co2_concentration)
    plt.grid(True)
    plt.xlabel("Year")
    plt.ylabel("Concentration [ppmv]")
    #plt.show()


#############################################################################
#For NOx emissions
#Assumption for s=1

t_test = np.arange(0,8974800,1)          #s - 1 year  - 4x500 missions
Type = 2
#1-NOx -> CH4
#2-NOx -> O3L
Configuration = 0
#0 - kerosene
#1 - Hydrogen Combustion

tau_i = 12*365*24*60    #years

if Type == 1:
    A_i = 5.16 * 10 ** -13  # W/m2/kg
if Type == 2:
    A_i = -1.21 * 10 ** -13  # W/m2/kg

if Configuration ==0:
    massflow=100            #kg/s
    EI_NOx = 0.014          #kg/kg
if Configuration==1:
    massflow = 8.5          # kg/s
    EI_NOx = 0.0028         # kg/kg
if Configuration==2:
    massflow = 186.7        # kg/s
    EI_NOx = 0.014          # kg/kg
if Configuration==3:
    massflow = 8.5          # kg/s
    EI_NOx = 0.014          # kg/kg
if Configuration==4:
    massflow = 8.5          # kg/s
    EI_NOx = 0.014          # kg/kg


def emissions_NOx(EI_NOx, massflow, t_year):
    output_NOx = EI_NOx*massflow*t_year
    return output_NOx

emissions_NOx(EI_NOx=EI_NOx, massflow=massflow, t_year=t_test)
#print("emissions NOx=", emissions_NOx(EI_NOx=EI_NOx, massflow=massflow, t_year=t_test))
emissionsNOx = np.sum(emissions_NOx(EI_NOx=EI_NOx, massflow=massflow, t_year=t_test))
print("Sum of NOx emissions=",emissionsNOx,"kg")


def Gfactor (t_year, A_i, tau_i):
    G = A_i*np.exp(-t_year/tau_i)
    return G

print("Gfactor=",Gfactor(t_year=t_test, A_i=A_i, tau_i=tau_i))
print(np.sum(Gfactor(t_year=t_test, A_i=A_i, tau_i=tau_i)))


def concentration_NOx(emissionsNOx, t_test, A_i,tau_i):
    """
    Estimates CO2 concentration in atmosphere in ppmv.

    :param t_year: Array with time, typically years
    :type t_year: np.ndarray
    :param emissions_co2: CO2 emissions per year in t_year, in kg
    :type emissions_co2: np.ndarray
    :param co2_0: Pre-industrial or background CO2 concentration in
        ppmv
    :type co2_0: float
    :param alpha_i: Decay factor in ppbv/Tg
    :param tau_i: Lifetime in years
    :return: CO2 concentration in atmosphere in ppmv
    """
    result = convolve(signal=emissionsNOx,response=Gfactor(t_year=t_test, A_i=A_i, tau_i=tau_i))
    return result

print("RFNOx =", np.sum(concentration_NOx(emissionsNOx=emissionsNOx, t_test=t_test, A_i=A_i, tau_i=tau_i)), "W/m^2")


'''
# Data CO2 from Sausen and Schumann (2000)
t_test = np.arange(1940, 2031, 1)
e_co2_test = []
fal_1940_1995 = np.array(
    [7.7, 8.3, 9.0, 9.7, 10.5, 11.3, 12.3, 13.2, 14.3, 15.4,
        16.7, 18.0, 19.4, 21.0, 22.7, 24.5, 26.5, 28.6, 30.9, 33.3,
         36.0, 39.5, 43.5, 45.9, 48.0, 51.3, 55.6, 65.6, 74.3, 77.8,
         78.0, 90.0, 96.0, 99.4, 96.0, 96.1, 96.4, 102.1, 105.7, 110.1,
         110.9, 109.3, 110.5, 112.0, 119.5, 123.4, 129.9, 135.6, 141.4, 146.5,
         146.9, 143.4, 142.0, 144.1, 150.0, 154.3])
for i, ti in enumerate(t_test):
    if ti <= 1995:
        e_co2_test.append(fal_1940_1995[i])
    elif 1995 < ti <= 2015:
        emi = ((278.7 - 154.3) / (2015 - 1995)) * (ti - 1995) + 154.3
        e_co2_test.append(emi)
    elif 2015 < ti <= 2050:
        emi = ((405.1 - 278.7) / (2050 - 2015)) * (ti - 2015) + 278.7
        e_co2_test.append(emi)
    elif 2050 < ti <= 2100:
        emi = ((666.2 - 405.1) / (2100 - 2050)) * (ti - 2050) + 405.1
        e_co2_test.append(emi)

e_co2_test = np.array(e_co2_test) * 1e9  # convert to kilograms

test_NOx_concentration = concentration_NOx(
        t_year=t_test,
        emissions_co2=e_co2_test)

def Convolution_RF_NOx():
        """
        Estimates NOx RF

        :param
        :type
        :param
        :type
        :param
        :type
        :param
        :param
        :return:
        """

delta = convolve(signal=,
                 response=)


def Gfactor(t_years,
                           alpha_i=np.array(
                               [0.067, 0.1135, 0.152, 0.097, 0.041]),
                           tau_i=np.array([np.inf, 318.8, 79.8, 18.8, 1.7])):
    """
    Delta CO2 concentration response due to CO2 emissions.

    :param t_years: Array with time, typically years
    :type t_years: np.ndarray
    :param alpha_i: Decay factor in ppbv/Tg
    :param tau_i: Lifetime in years
    :return: Delta CO2 concentration response function
    """


    tn = t_years - t_years[0]  # + 1
    alpha_in = alpha_i * 1e-3 * 1e-9  # convert from ppbv/Tg to ppmv/kg
    try:
        one = np.ones((5, tn.size)) * -1.
    except AttributeError:
        one = np.ones((5, 1)) * -1.
    one = one * tn
    one = one.transpose()
    one = one / tau_i
    one = np.exp(one)

    one = one * alpha_in

    one = np.sum(one, axis=1)

    return one






'''

'''
#1-NOx long CH4
#2-NOx long O3L
#3-CO2
Type = 2

#NOx
tau_CH = 12         #years
t = 13              #years


if Type==1:
    A = -5.16*10**-13
    E = 0.119*t
    RF = (A * np.exp(-t / tau_CH)) * (np.exp(t / tau_CH) * E - (tau_CH ** 2) * (np.exp(t / tau_CH) - np.exp(0)))
    print("RF_NOXlong_CH4=", RF)
if Type==2:
    A = -1.21*10**-13
    E = 0.119*t
    RF = (A * np.exp(-t / tau_CH)) * (np.exp(t / tau_CH) * E - (tau_CH ** 2) * (np.exp(t / tau_CH) - np.exp(0)))
    print("RF_NOXlong_O3L=", RF)
if Type==3:
'''
