import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ja')


def sigma_real(s, ge, gm, gt, mp, mn):
    r"""
    In the case of real dark matter, computes the cross section for the scattering of a :math:`\mu`-neutrino off a cosmic neutrino background neutrino into two scalar dark matter particles.

    Parameters
    ----------
    s : double
        centre of mass energy [:math:`\mathrm{MeV}^2`]
    ge, gm, gt : double
        coupling between electron/mu/tau neutrino and dark sector
    mp, mn : double, double
        masses of the dark scalar and dark fermion respectively [MeV]

    Returns
    -------
    sigma : double
        the resulting cross section [:math:`\mathrm{MeV}^{-2}`]

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from crossSection import sigma_real; import numpy as np;

       In [2]: sigma_real(s=6.0, ge=0, gm=np.power(10.0, -2), gt=3*np.power(10.0, -1), mp=1.0, mn=2.0)
    """
    try:
        prefactor = 2*np.power(gm, 2)*sum(np.power([ge, gm, gt], 2))*np.power(mn, 2)/(32.0*np.pi)
        first_term = 2*np.sqrt((0.25*s - np.power(mp, 2))/s)/(np.power(np.power(mp, 2) - np.power(mn, 2) - 0.5*s, 2) - s*(0.25*s - np.power(mp, 2)))
        second_term = np.log((np.power(mp, 2) - np.power(mn, 2) - 0.5*s - np.sqrt(s*(0.25*s - np.power(mp, 2))))/(np.power(mp, 2) - np.power(mn, 2) - 0.5*s + np.sqrt(s*(0.25*s - np.power(mp, 2)))))/(s*(np.power(mp, 2) - np.power(mn, 2) - 0.5*s))
        sigma = prefactor*(first_term + second_term)
        return sigma
    except:
        sigma = np.nan
        return sigma


def sigma_without_prefactor(s, mp, mn):
    r"""
    Parameters
    ----------
    s : double
        centre of mass energy [:math:`\mathrm{MeV}^2`]
    ge, gm, gt : double
        coupling between electron/mu/tau neutrino and dark sector
    mp, mn : double, double
        masses of the dark scalar and dark fermion respectively [MeV]

    Returns
    -------
    sigma : double
        the resulting cross section [:math:`\mathrm{MeV}^{-2}`]

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from crossSection import sigma_without_prefactor

       In [2]: sigma_without_prefactor(s=6.0, mp=1.0, mn=2.0)
    """
    try:
        first_term = 2*np.sqrt((0.25*s - np.power(mp, 2))/s)/(np.power(np.power(mp, 2) - np.power(mn, 2) - 0.5*s, 2) - s*(0.25*s - np.power(mp, 2)))
        second_term = np.log((np.power(mp, 2) - np.power(mn, 2) - 0.5*s - np.sqrt(s*(0.25*s - np.power(mp, 2))))/(np.power(mp, 2) - np.power(mn, 2) - 0.5*s + np.sqrt(s*(0.25*s - np.power(mp, 2)))))/(s*(np.power(mp, 2) - np.power(mn, 2) - 0.5*s))
        sigma = first_term + second_term
        return sigma
    except:
        sigma = 0
        return sigma

def PMNS_matrix(t12, t23, t13, dcp):
    r"""
    Returns the absolute values of the 3x3 PMNS matrix as an array with rows (e, :math:`\mu`, :math:`\tau`) and columns (1, 2, 3)

    Returns 
    -------
    pmns_abs : np.array
        PMNS matrix, see `here <https://en.wikipedia.org/wiki/Pontecorvo%E2%80%93Maki%E2%80%93Nakagawa%E2%80%93Sakata_matrix>`_
    
    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from crossSection import PMNS_matrix

       In [2]: PMNS_matrix(t12 = 33.63, t23 = 47.2, t13 = 8.54, dcp = 234)
    """
    # Define the various sines, cosines of mixing angles
    c12 = np.cos(np.deg2rad(t12))
    c13 = np.cos(np.deg2rad(t13))
    c23 = np.cos(np.deg2rad(t23))
    s12 = np.sin(np.deg2rad(t12))
    s13 = np.sin(np.deg2rad(t13))
    s23 = np.sin(np.deg2rad(t23))

    # Initialise empty complex array
    pmns_abs = np.empty([3,3])

    # Fill the array
    pmns_abs[0][0] = np.absolute(c12*c13)
    pmns_abs[1][0] = np.absolute(-s12*c23 - c12*s23*s13*np.exp(dcp*1.0j))
    pmns_abs[2][0] = np.absolute(s12*s23 - c12*c23*s13*np.exp(dcp*1.0j))
    pmns_abs[0][1] = np.absolute(s12*c13)
    pmns_abs[1][1] = np.absolute(c12*c23 - s12*s23*s13*np.exp(dcp*1.0j))
    pmns_abs[2][1] = np.absolute(-c12*s23 - s12*c23*s13*np.exp(dcp*1.0j))
    pmns_abs[0][2] = np.absolute(s13*np.exp(-dcp*1.0j))
    pmns_abs[1][2] = np.absolute(s23*c13)
    pmns_abs[2][2] = np.absolute(c23*c13)

    return pmns_abs

def mfp_gpc(n, sigma, cm):
    r"""
    Returns the mean free path in Gpc. If the cross-section vanishes, the division will raise an exception and a large value is returned.

    Parameters
    ----------
    n     : double
    	number density of target in :math:`\mathrm{cm}^{-3}`
    sigma : double
    	cross-section in :math:`\mathrm{cm}^{2}`
    cm    : double
    	conversion factor from cm to Gpc

    Returns
    -------
    mfp_gpc : double
    	mean free path in Gpc

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from crossSection import mfp_gpc

       In [2]: mfp_gpc(n=1, sigma=np.power(10.0, -32), cm=np.power(10.0, 21))
    """
    try:
        # 1.903 factor accounts for the redshift up to z = 0.3365
        mfp_gpc = (1.0/1.903)*np.power(n*sigma, -1)*cm
        return mfp_gpc
    except:
        # Likely that sigma vanishes as centre of mass energy below threshold
        # Return standard model value of 10^{11} Gpc
        return np.power(10.0, 11)

def neutrino_masses(m_min, hierarchy):
    r"""
    Takes the lightest mass and the hierarchy as parameters, returns an array [m1, m2, m3] with the neutrino masses in each of the two heirarchies.

    Parameters
    ----------
    m_min : double
    	mass of the lightest neutrino species in eV
    hierarchy : string
    	choice of hierarchy, either 'normal' or 'inverted' (case insensitive)

    Returns
    -------
    m_array : np.array
    	masses of the mass eigenstates [m1, m2, m3]

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from crossSection import neutrino_masses

       In [2]: neutrino_masses(m_min=0.03, hierarchy='normal')
    """
    if hierarchy.lower() not in ['normal', 'inverted']:
        print('WARNING: Choice of hierarchy not recognised, please choose either \'normal\' or \'inverted\', defaulting to NORMAL hierarchy\n')
        hierarchy = 'normal'

    if hierarchy.lower() == 'normal':
        print('Computing neutrino masses in NORMAL hierarchy')
        m1 = m_min
        m2 = np.sqrt(7.37*np.power(10.0, -5) + np.power(m1, 2))
        m3 = np.sqrt(2.50*np.power(10.0, -3) + 0.5*(np.power(m1, 2) + np.power(m2, 2)))
        mass_sum = m1 + m2 + m3
        print('-----------------------')
        print('m1           = {0:.3f} eV'.format(m1))
        print('m2           = {0:.3f} eV'.format(m2))
        print('m3           = {0:.3f} eV'.format(m3))
        print('m1 + m2 + m3 = {0:.3f} eV'.format(mass_sum))
        print('-----------------------\n')

    elif hierarchy.lower() == 'inverted':
        print('Computing neutrino masses in INVERTED hierarchy')
        m3 = m_min
        m2 = np.sqrt(0.5*7.37*np.power(10.0, -5) + np.power(m3, 2) + 2.46*np.power(10.0, -3))
        m1 = np.sqrt(np.power(m2, 2) - 7.37*np.power(10.0, -5))
        mass_sum = m1 + m2 + m3
        print('-----------------------')
        print('m1           = {0:.3f} eV'.format(m1))
        print('m2           = {0:.3f} eV'.format(m2))
        print('m3           = {0:.3f} eV'.format(m3))
        print('m1 + m2 + m3 = {0:.3f} eV'.format(mass_sum))
        print('-----------------------\n')

    return np.array([m1, m2, m3])


def sigma_with_masses(s, ge, gm, gt, mp, mn, pmns):
    r"""
    Takes an np.array, s, as input which contains the centre of mass energies for the three neutrino mass eigenstates. Also, given the flavour couplings, computes the mass couplings gi. Outputs the sum across the three mass eigenstates.

    Parameters
    ----------
    s : np.array
    	np.array of length 3 with centre of mass energies
    ge, gm, gt  : double
    	coupling between electron/mu/tau neutrino and dark sector
    mp, mn : double
    	masses of the dark scalar and dark fermion respectively
    pmns : np.arrray
    	absolute value of the PMNS matrix to compute mass couplings

    Returns
    -------
    sigma_total : float
    	returns sum of the cross sections evaluated at the different centre of mass energies

    Examples
    --------
    .. ipython::
       :okwarning:

       @suppress
       In [1]: from crossSection import sigma_with_masses, PMNS_matrix

       In [2]: sigma_with_masses(s=np.array([10**2, 8**2, 12**2]), ge=1, gm=1, gt=1, mp=1, mn=2, pmns=PMNS_matrix(t12 = 33.63, t23 = 47.2, t13 = 8.54, dcp = 234))
    """
    if len(s) != 3:
        print('ERROR: len(s) is not equal to 3')
    # Compute mass couplings
    g1 = pmns[0][0]*ge + pmns[1][0]*gm + pmns[2][0]*gt
    g2 = pmns[0][1]*ge + pmns[1][1]*gm + pmns[2][1]*gt
    g3 = pmns[0][2]*ge + pmns[1][2]*gm + pmns[2][2]*gt

    # Compute prefactors (factor of 2 for antiparticles, complex factor of 3 in main code)
    prefactor1 = 2*np.power(gm, 2)*np.power(g1, 2)*np.power(mn, 2)/(32.0*np.pi)
    prefactor2 = 2*np.power(gm, 2)*np.power(g1, 2)*np.power(mn, 2)/(32.0*np.pi)
    prefactor3 = 2*np.power(gm, 2)*np.power(g1, 2)*np.power(mn, 2)/(32.0*np.pi)

    # Compute contribution from each mass state
    sigma1 = prefactor1*sigma_without_prefactor(s[0], mp, mn)
    sigma2 = prefactor2*sigma_without_prefactor(s[1], mp, mn)
    sigma3 = prefactor2*sigma_without_prefactor(s[2], mp, mn)

    return sigma1 + sigma2 + sigma3

if __name__ == '__main__':
    print('Nothing to compute')
    # # Constants and Units
    # eV  = np.power(10.0, -6)
    # MeV = np.power(10.0, 6)*eV
    # GeV = np.power(10.0, 9)*eV
    # TeV = np.power(10.0, 12)*eV
    # PeV = np.power(10.0, 15)*eV
    # ge  = np.power(10.0, -3)
    # gm  = np.power(10.0, -3)
    # gt  = np.power(10.0, -3)
    #
    # mp  = 1*MeV
    # mn  = 10*MeV
    #
    # # Arrays
    # energy_arr = np.linspace(mp, 10*mn, 500)     # Centre of Mass energies
    # s_arr = 4*np.power(energy_arr, 2)
    # sigma_list = []
    # for s in s_arr:
    #     sigma_list.append(sigma_real(s, ge, gm, gt, mp, mn))
    # sigma_arr = np.array(sigma_list)
    #
    # # Plotting
    # plt.figure()
    # plt.plot(energy_arr, sigma_arr, c='r')
    # axes = plt.axis()
    # new_axes = [0, max(energy_arr), 0, axes[3]]
    # plt.axis(new_axes)
    # ax = plt.gca()
    # ax.axvline(mp, linestyle='-.', color='k')
    # plt.xlabel(r'Centre of Mass Energy, $E / \mathrm{MeV}$')
    # plt.ylabel(r'Cross Section, $\sigma / \mathrm{MeV}^{-2}$')
    # style = dict(size=6, color='r')
    # plt.text(energy_arr[300], 0.3*sigma_arr[400], r'$g_{e, \mu, \tau} = 10^{-3}, m_\phi = 5\,\mathrm{MeV}, m_N = 10\,\mathrm{MeV}$', **style)
    # plt.savefig('exampleSigma.pdf')
