from crossSection import sigma_with_masses, PMNS_matrix, neutrino_masses, mfp_gpc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import random
import warnings
import sys
import os
warnings.filterwarnings("ignore")
plt.style.use('ja')
data_dir = '/Users/james/allMyStuff/Neutrinos/Constraints/highenergy/'

if __name__ == '__main__':
    print('-----------------------------------')
    print('     Running include_masses.py     ')
    print('-----------------------------------\n')
    # Define blazar neutrino energy
    Enu_TeV = np.array([2000.])

    params = {'mp': float(sys.argv[1]), 'min_nu_mass': float(sys.argv[2]), 'hierarchy': sys.argv[3]}

    # Define minimum neutrino mass and calculate others in hierarchy, returns array in eV
    min_nu_mass = params['min_nu_mass'] # 0.03 # eV
    hierarchy = params['hierarchy'] # 'normal'
    nu_masses = neutrino_masses(min_nu_mass, hierarchy)

    # Calculate c.o.m energies as np.array
    Ecom_MeV = np.sqrt(0.5*nu_masses*Enu_TeV)
    s_arr = 4 * np.power(Ecom_MeV, 2)
    print('s =', np.sqrt(s_arr))

    print('-----------------------------------')
    print('            Parameters             ')
    print('-----------------------------------\n')

    # Fix ge and gt to be maximal in the complex case
    ge = 0 # 3*np.power(10.0, -3)
    gt = 3*np.power(10.0, -1) # 3*np.power(10.0, -1)
    #print('ge = {0:.5f}, gt = {0:.5f}'.format(ge, gt))
    print('ge = {0:.4f}'.format(ge))
    print('gt = {0:.4f}'.format(gt))

    # Fix mp for now, TODO: Vary later
    mp = params['mp'] # 1.0 # MeV
    print('mscalar = {} MeV'.format(mp))

    # Fix neutrino number density and mass (multiply by 3 in complex case assuming small splitting)
    n_nu = 340 # cm^-3
    n_eff = (340/6.0) * 3.0

    # Conversion factors
    cm = 3.240755 * np.power(10.0, -28) # Gpc
    MeV = 8065.54429 * np.power(10.0, 6) # cm^-1

    # Distance to blazar
    D_blazar = 1.3 # Gpc
    print('Distance to blazar: {} Gpc\n'.format(D_blazar))

    # PMNS matrix

    t12 = 33.63
    t23 = 47.2
    t13 = 8.54
    dcp = 234
    pmns = PMNS_matrix(t12, t23, t13, dcp)
    print('PMNS Matrix')
    print(pmns)
    print('\n')

    # Limits

    mn_max = 20.0 # MeV
    g_min = -5.0
    g_max = -1.0
    axis_min = -4.0
    axis_max = -1.0

    # Experiment loops
    ng = 200
    nm = 200
    count = 1
    gmin_arr = np.array([])
    mn_arr = np.array([])

    print('Starting to search parameter space...\n')
    for mn in np.linspace(0.0, mn_max, nm):
        for gm in np.logspace(g_min, g_max, ng):
            print('Searched (' + str(count) + '/' + str(ng*nm) + ')', end='\r')
            sigma_cm = sigma_with_masses(s_arr, ge, gm, gt, mp, mn, pmns)*np.power(MeV, -2)
            mfp = mfp_gpc(n_eff, sigma_cm, cm)

            if mfp < D_blazar:
                gmin_arr = np.append(gmin_arr, gm)
                mn_arr = np.append(mn_arr, mn)
                break

            count += 1

    save = False
    if save:
        filename = 'mp' + str(mp) + 'numin' + str(min_nu_mass) + hierarchy.upper() + '.csv'
        data = pd.DataFrame({'gmin_arr': gmin_arr, 'mn_arr': mn_arr})
        data.to_csv(data_dir + filename)
        print('--- Saved to ' + data_dir  + filename + ' ---')

    # Constraint Plotting
    else:
        plt.figure()
        plt.rcParams.update({'font.size': 22})

        plt.plot([mp, mp], [np.power(10.0, -4), np.power(10.0, -2)], c='k', linewidth=1.0)
        plt.plot([mp, 20.0], [np.power(10.0, -2), np.power(10.0, -2)], c='k', linewidth=1.0)
        #plt.plot([10.0, 10.0], [np.power(10.0, -2), np.power(10.0, -4)], c='k', linewidth=1.0)

        plt.loglog(mn_arr, gmin_arr, linewidth=1.0, linestyle='-', marker='', markersize=0.0, markerfacecolor='r', alpha = 1.0, markeredgewidth=0.0)

        upper_limit = np.empty(len(gmin_arr))
        upper_limit.fill(np.power(10.0, g_max))

        plt.fill_between(mn_arr, gmin_arr, upper_limit, alpha=0.2, edgecolor='r', facecolor='k', linewidth=0.0)

        #plt.fill_between([10.0, mn_max], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

        plt.fill_between([0.0, mp], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)
        # plt.fill_betweenx([np.power(10.0, -5), 3*np.power(10.0,-4)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

        plt.fill_betweenx([np.power(10.0, -2), np.power(10.0, g_max)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

        axes = plt.axis()
        plt.axis([1.0, mn_max, np.power(10.0, axis_min), np.power(10.0, axis_max)])

        plt.xlabel(r'$m_N / \mathrm{MeV}$')
        plt.ylabel(r'$g_\mu$')
        plt.savefig('constraintsExample.pdf')
        print('\n\n--- Saved constraintsExample.pdf ---')
        os.system('open constraintsExample.pdf')
