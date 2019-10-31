import numpy as np
import pandas as pd


def integrate_trap(x, y):
    r"""
    Rudimentary application of the trapezium rule for discrete data sets x, y

    Parameters
    ----------
    x, y : np.array or list
    	(x, y) data to numerically integrate, applies a standard trapezium rule.

    Returns
    -------
    sum : float
    	the result of the numerical integration

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from integration import integrate_trap

       In [2]: integrate_trap(x=np.linspace(0, 1), y=np.power(np.linspace(0, 1), 2))
    """
    n = len(x)
    sum = 0
    for k in range(1, n):
        sum += 0.5*(y[k - 1] + y[k])*(x[k] - x[k - 1])
    return sum


def luminosity(d, flux):
    r"""
    Given a flux, returns the luminosity using :math:`L = 4\pi D^2 F`

    Parameters
    ----------
    d    : float
    	distance to the source
    flux : float
    	flux of particles etc. coming from the source

    Returns
    -------
    lum : float
    	luminosity

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from integration import luminosity

       In [2]: luminosity(d=1.0, flux=1.0)
    """
    return 4*np.pi*np.power(d, 2)*flux


# Global definitions

if __name__ == '__main__':
    Gpc = 3.085678*np.power(10.0, 27)
    TeV = np.power(10.0, -12)
    dL = 1.3*Gpc
    data_dir = '/Users/james/allMyStuff/Neutrinos/IceCubeData/luminosity/'

    print('\nHAWC Data\n')
    hawc_data = pd.read_csv(data_dir+'hawc.csv', header=0)
    hawc_E = np.array(hawc_data['E'])
    hawc_vf = np.array(hawc_data['vf'])
    print('- Computing luminosity from HAWC data in range {} TeV to {} TeV \n'.format(round(hawc_E[0]*TeV, 2), round(hawc_E[-1]*TeV), 2))
    hawc_logE = np.log(hawc_E)

    hawc_flux = integrate_trap(hawc_logE, hawc_vf)
    print('- Flux from HAWC data: {:0.2e} erg/s/cm\u00b2 \n'.format(hawc_flux))

    hawc_luminosity = luminosity(dL, hawc_flux)
    print('- Luminosity from HAWC data: {:0.2e} erg/s \n'.format(hawc_luminosity))


    print('\nNeutrino Data\n')
    nu_data = pd.read_csv(data_dir+'nu.csv', header=0)
    nu_E = np.array(nu_data['E'])
    nu_vf_max = np.array(nu_data['vf_max'])
    nu_vf = np.array(nu_data['vf'])
    nu_vf_min = np.array(nu_data['vf_min'])
    print('- Computing luminosity from Neutrino data in range {} TeV to {} TeV \n'.format(round(nu_E[0]*TeV, 2), round(nu_E[-1]*TeV), 2))
    nu_logE = np.log(nu_E)

    nu_flux_max = integrate_trap(nu_logE, nu_vf_max)
    nu_flux = integrate_trap(nu_logE, nu_vf)
    nu_flux_min = integrate_trap(nu_logE, nu_vf_min)
    print('- Avg. flux from Neutrino data: {:0.2e} erg/s/cm\u00b2 \n'.format(nu_flux))

    nu_luminosity_max = luminosity(dL, nu_flux_max)
    nu_luminosity = luminosity(dL, nu_flux)
    nu_luminosity_min = luminosity(dL, nu_flux_min)
    print('- Mean luminosity from Neutrino data: {:0.2e} erg/s \n'.format(nu_luminosity))
    print('- Min luminosity from Neutrino data: {:0.2e} erg/s \n'.format(nu_luminosity_min))
    print('- Max luminosity from Neutrino data: {:0.2e} erg/s \n'.format(nu_luminosity_max))
