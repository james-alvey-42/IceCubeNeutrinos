if __name__ == '__main__':
    from crossSection import sigma_real
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import sys
    import random
    import warnings
    warnings.filterwarnings("ignore")
    plt.style.use('ja')
    from numpy.random import normal, sample
    Enu_TeV = np.array([290.])
    ge = 3*np.power(10.0, -3)
    gt = 3*np.power(10.0, -1)
    gm = np.power(10.0, -2)
    mp = 0.5
    mp_MeV = mp*np.power(10.0, 6)
    nu_mass = 0.15
    mn = 5

    Ecom_nu_MeV = np.sqrt(0.5*nu_mass*Enu_TeV)
    Ecom_phi_MeV = np.sqrt(0.5*mp_MeV*Enu_TeV)

    s_nu_sample = 4 * np.power(Ecom_nu_MeV, 2)
    s_phi_sample = 4 * np.power(Ecom_phi_MeV, 2)
    sigma_nu = sigma_real(s_nu_sample, ge, gm, gt, mp, mn)
    sigma_phi = sigma_real(s_phi_sample, ge, gm, gt, mp, mn)

    print('Scattering off Neutrino:', sigma_nu[0])
    print('Scattering off Scalar:', sigma_phi[0])
