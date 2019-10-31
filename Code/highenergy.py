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
	Enu_TeV_arr = [290., 3000., 6000.]
	labels = ['290 TeV', '3 PeV', '6 PeV']
	colors = ['#357DED', '#0D0221', '#0D0221']
	lss = ['-', '--', ':']
	plt.figure(figsize=(8, 5))
	for idx, Enu_TeV in enumerate(Enu_TeV_arr):
		min_nu_mass = 0.03 # GeV
		hierarchy = 'normal'
		nu_masses = neutrino_masses(min_nu_mass, hierarchy)
		Ecom_MeV = np.sqrt(0.5*nu_masses*Enu_TeV)
		s_arr = 4 * np.power(Ecom_MeV, 2)
		ge = 0
		gt = 3*np.power(10.0, -1)
		gm = np.power(10.0, -2)
		# Fix neutrino number density and mass (multiply by 3 in complex case assuming small splitting)
		n_nu = 340 # cm^-3
		n_eff = (340/6.0) * 3.0

		# Conversion factors
		cm = 3.240755 * np.power(10.0, -28) # Gpc
		MeV = 8065.54429 * np.power(10.0, 6) # cm^-1

		# Distance to blazar
		D_blazar = 1.3 # Gpc

		# PMNS matrix

		t12 = 33.63
		t23 = 47.2
		t13 = 8.54
		dcp = 234
		pmns = PMNS_matrix(t12, t23, t13, dcp)

		mn = np.linspace(0.0, 15.0, 500)
		mp = np.linspace(0.0, 15.0, 500)

		MN, MP = np.meshgrid(mn, mp)

		sigma_cm = sigma_with_masses(s_arr, ge, gm, gt, MP, MN, pmns)*np.power(MeV, -2)
		mfp = mfp_gpc(n_eff, sigma_cm, cm)
		region = (MN <= MP)
		mfp[region] = np.ma.masked

		ctr = plt.contour(MP, MN, mfp, 
			colors=colors[idx], 
			levels=[D_blazar],
			linewidths=0.0)
		mp_trace, mn_trace = ctr.allsegs[0][0].T
		mask = (mp_trace < 0.9*mn_trace)
		plt.plot(mp_trace[mask], mn_trace[mask], 
			c=colors[idx],
			ls=lss[idx],
			label=labels[idx])
		if idx == 0:
			plt.plot([mp_trace[0], mp_trace[0]], [mn_trace[0], mp_trace[0]],
				c=colors[idx],
				ls=lss[idx],)
		plt.fill(np.append(mp_trace, [0, mp_trace[0]]), np.append(mn_trace, [0, mp_trace[0]]), 
			color=colors[idx], 
			alpha=0.1)
	plt.plot(np.linspace(0.1, 10), np.linspace(0.1, 10),
		color='k',
		ls='-',
		lw=0.3)
	plt.text(np.power(10.0, -0.7), np.power(10.0, -0.6), r'$m_N > m_\delta$', rotation=27.0)
	plt.xlabel(r'$m_\delta \, \mathrm{[MeV]}$')
	plt.ylabel(r'$m_N \, \mathrm{[MeV]}$')
	plt.plot([3.9, 3.9], [0.0, 30.0], 
		label='',
		c='k',
		ls='-')
	plt.xscale('log')
	plt.yscale('log')
	plt.text(4.3, np.power(10.0, -0.8), r'$\mathrm{N}_{\mathrm{eff}}$', rotation=90.0)
	plt.plot([6.74, 6.74], [0.0, 30.0], 
		label='',
		c='k',
		ls='-')
	plt.text(7.4, np.power(10.0, 0.5), r'$\mathrm{BBN} + \mathrm{Planck} + \mathrm{N}_{\mathrm{eff}} + Y_p$', rotation=90.0)
	plt.text(np.power(10.0, -0.9), np.power(10.0, 1.1), r'$g_\mu = 10^{-2}, m_\nu^{\mathrm{min}} = 0.03 \, \mathrm{eV}$')
	plt.text(np.power(10.0, -0.9), np.power(10.0, 1.3), r'$\mathrm{Normal}\,\,\mathrm{Hierarchy}$')
	#plt.plot(np.linspace(0, 12), np.linspace(0, 12))
	plt.legend(loc='lower center', fontsize=14, title='Neutrino Energy', title_fontsize=10)
	plt.xlim([np.power(10.0, -1), np.power(10.0, 1)])
	plt.ylim([np.power(10.0, -1), 30.0])
	ax = plt.gca()
	ax.tick_params(which='minor', length=1.5)
	plt.savefig('mpmnconstraints_test.pdf')
	
