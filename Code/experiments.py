from crossSection import sigma_real
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import random
import warnings
warnings.filterwarnings("ignore")
plt.style.use('ja')
data_dir = '/Users/james/allMyStuff/Neutrinos/IceCubeData/'

if __name__ == '__main__':
    from numpy.random import normal, sample
    number_of_events = int(np.floor(normal(13.0, 5.0, 1)[0]))
    print('Number of Events: ' + str(number_of_events) + '\n')

    txs_df = pd.read_csv(filepath_or_buffer=data_dir+'txsEvents.csv', header=0, sep='   ', engine='python')
    indices = list(range(0, len(txs_df)))
    mjd_arr = np.array(txs_df['MJD'])
    logE_arr = np.array(txs_df['log10(Ereco)'])

    # Take random sample of #(number_of_events) TXS Events
    indices_sample = random.sample(indices, number_of_events)
    mjd_sample = mjd_arr[indices_sample]
    logE_sample = logE_arr[indices_sample]

    # Data Bar Chart to show sample
    plt.figure()
    plt.bar(mjd_arr, logE_arr,
            facecolor='k',
            width=1.0,
            linewidth=0,
            align='edge')
    plt.bar(mjd_sample, logE_sample,
            facecolor='r',
            width=1.0,
            linewidth=0,
            align='edge')
    # plt.plot(mjd_arr, logE_arr, 'r-.', linewidth=0.5)
    plt.xlabel('Modified Julian Date')
    plt.ylabel(r'$\log_{10} \left[ \frac{\hat{E}_{\mu}}{\textrm{\tiny GeV}} \right]$')
    plt.savefig('sample.pdf')
    print('--- Saved sample.pdf ---\n')

    Emu_TeV = 10**(-3) * 10**logE_sample
    # Enu_TeV = 1.92 * Emu_TeV**(1.14)    # From Kelly et al. (2018)
    Enu_TeV = np.array([290.])

    print('--- NOTE: Using PeV Neutrino ---\n')

    nu_mass = 0.01 # eV

    # Note that TeV * eV = (MeV)^2
    Ecom_MeV = np.sqrt(0.5*nu_mass*Enu_TeV)
    s_sample = 4 * np.power(Ecom_MeV, 2)

    # Fix ge and gt to be maximal in the complex case
    ge = 3*np.power(10.0, -3)
    gt = 3*np.power(10.0, -1)

    # Fix mp for now
    # TODO: Vary later
    mp = 0.5 # MeV

    # Fix neutrino number density and mass
    n_nu = 340 # cm^-3
    n_eff = 340/6.0 * 3.0

    # Conversion factors
    cm = 3.240755 * np.power(10.0, -28) # Gpc
    MeV = 8065.54429 * np.power(10.0, 6) # cm^-1

    # Distance to blazar
    D_blazar = np.power(5.0, 0) # Gpc

    # Limits

    mn_max = 12.0 # MeV
    g_min = -5.0
    g_max = 1.0
    axis_min = -4.0
    axis_max = -1.0

    # Experiment loops
    ng = 200
    nm = 200
    count = 1
    gmin_arr = np.array([])
    mn_arr = np.array([])
    for mn in np.linspace(0.0, mn_max, nm):
        for gm in np.logspace(g_min, g_max, ng):
            print('Searched (' + str(count) + '/' + str(ng*nm) + ')', end='\r')
            sigma_sample_cm = sigma_real(s_sample, ge, gm, gt, mp, mn)*np.power(MeV, -2)
            mfp_gpc = (1.0/1.903)*np.power(n_eff*sigma_sample_cm, -1)*cm
            if not all(l > D_blazar for l in mfp_gpc):
                if not np.isnan(mfp_gpc).all():
                    gmin_arr = np.append(gmin_arr, gm)
                    mn_arr = np.append(mn_arr, mn)
                    break
            count += 1

    save = False
    if save:
        filename = 'mp' + str(mp) + 'mnu' + str(nu_mass) + '.csv'
        data = pd.DataFrame({'gmin_arr': gmin_arr, 'mn_arr': mn_arr})
        data.to_csv(data_dir + 'complex/' + filename)
        print('--- Saved to ' + data_dir + 'complex/' + filename + ' ---')

    # Constraint Plotting
    else:
        plt.figure()

        plt.semilogy(mn_arr, gmin_arr, linewidth=0.0, marker='', markersize=0.0, markerfacecolor='r', alpha = 1.0, markeredgewidth=0.0)

        upper_limit = np.empty(len(gmin_arr))
        upper_limit.fill(np.power(10.0, g_max))

        plt.fill_between(mn_arr, gmin_arr, upper_limit, alpha=0.4, edgecolor='r', facecolor='r', linewidth=0.0)

        plt.fill_between([10.0, mn_max], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

        plt.fill_between([0.0, mp], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)
        # plt.fill_betweenx([np.power(10.0, -5), 3*np.power(10.0,-4)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

        plt.fill_betweenx([np.power(10.0, -2), np.power(10.0, g_max)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

        axes = plt.axis()
        plt.axis([0.0, mn_max, np.power(10.0, axis_min), np.power(10.0, axis_max)])

        plt.xlabel(r'$m_N / \mathrm{MeV}$')
        plt.ylabel(r'$g_\mu$')
        plt.savefig('constraintsExample.pdf')
        print('\n\n--- Saved constraints.pdf ---')
