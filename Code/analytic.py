import matplotlib.pyplot as plt
import numpy as np
from math import pi, log

def alpha(s, m, M):
    return m**2 - M**2 - 0.5*s

def beta(s, m, M):
    return (s*(s/4 - m**2))**(0.5)

def sigma(s, m, M, g):
    prefactor = ((s - 4*m**2)/s)**(0.5)*(g**4 * M**2)/(128*pi)
    term1 = prefactor*(2*beta(s, m, M)/(alpha(s, m, M)**2 - beta(s, m, M)**2))
    term2 = prefactor*(2/(alpha(s, m, M)*beta(s, m, M))*log((alpha(s, m, M) + beta(s, m, M))/(alpha(s, m, M) - beta(s, m, M)))
    return term1 + term2

if __name__ == '__main__':
    MeV = 10**6
    GeV = 10**9
    PeV = 10**15
    m = 1 MeV
    M = 10 MeV
    g = 10**(-3)
    root_s_array = np.linspace(0, 2*M, 1000)
    succesful_root_s_array = []
    cross_section_array = []
    for roots in root_s_array:
        try:
            temp = sigma(roots**2, m, M, g)
            cross_section_array.append(temp)
            succesful_root_s_array.append(roots)
        except:
            continue
    cross_section_array = np.array(cross_section_array)
    plt.plot(successful_root_s_array, cross_section_array)
    plt.show()
