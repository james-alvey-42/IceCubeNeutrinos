if __name__ == '__main__':
	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	import warnings
	warnings.filterwarnings("ignore")
	plt.style.use('ja')
	data_dir = '../IceCubeData/'

	mp = 1.0
	nu_mass = 0.15
	filename = 'mp' + str(mp) + 'mnu' + str(nu_mass) + '.csv'
	data = pd.read_csv(data_dir + 'complex/' + filename, index_col=0)
	print('--- Opened ' + data_dir + 'complex/' + filename + ' ---')
	gmin_arr1 = data['gmin_arr']
	mn_arr1 = data['mn_arr']

	mp = 0.5
	nu_mass = 0.1
	filename = 'mp' + str(mp) + 'mnu' + str(nu_mass) + '.csv'
	data = pd.read_csv(data_dir + 'complex/' + filename, index_col=0)
	print('--- Opened ' + data_dir + 'complex/' + filename + ' ---')
	gmin_arr2 = data['gmin_arr']
	mn_arr2 = data['mn_arr']

	mp = 0.5
	nu_mass = 0.05
	filename = 'mp' + str(mp) + 'mnu' + str(nu_mass) + '.csv'
	data = pd.read_csv(data_dir + 'complex/' + filename, index_col=0)
	print('--- Opened ' + data_dir + 'complex/' + filename + ' ---')
	gmin_arr3 = data['gmin_arr']
	mn_arr3 = data['mn_arr']

	mn_max = 12.0 # MeV
	g_min = -5.0
	g_max = -1.0
	axis_min = -4.0
	axis_max = -1.0

	plt.figure()

	plt.rcParams.update({'font.size': 22})

	plt.plot([mp, mp], [np.power(10.0, -4), np.power(10.0, -2)], c='k', linewidth=1.0)
	plt.plot([mp, 10.0], [np.power(10.0, -2), np.power(10.0, -2)], c='k', linewidth=1.0)
	plt.plot([10.0, 10.0], [np.power(10.0, -2), np.power(10.0, -4)], c='k', linewidth=1.0)

	plt.semilogy(mn_arr1, gmin_arr1, linewidth=1.0, linestyle='-', marker='', markersize=0.0, markerfacecolor='r', alpha = 1.0, markeredgewidth=0.0)

	upper_limit1 = np.empty(len(gmin_arr1))
	upper_limit1.fill(np.power(10.0, g_max))

	plt.fill_between(mn_arr1, gmin_arr1, upper_limit1, alpha=0.2, edgecolor='r', facecolor='k', linewidth=2.0)

	style = dict(size=15, color='r')
	plt.text(6.0, np.power(10.0, -2.6), r'$m_\nu = 0.15 \, \mathrm{eV}$', **style)

	plt.semilogy(mn_arr2, gmin_arr2, linewidth=1.0, linestyle='-', marker='', markersize=0.0, markerfacecolor='b', alpha = 1.0, markeredgewidth=0.0)

	upper_limit2 = np.empty(len(gmin_arr2))
	upper_limit2.fill(np.power(10.0, g_max))

	plt.fill_between(mn_arr2, gmin_arr2, upper_limit2, alpha=0.2, edgecolor='b', facecolor='k', linewidth=0.0)

	style = dict(size=15, color='b')
	plt.text(6.0, np.power(10.0, -2.9), r'$m_\nu = 0.10 \, \mathrm{eV}$', **style)

	plt.semilogy(mn_arr3, gmin_arr3, linewidth=1.0, linestyle='-', marker='', markersize=0.0, markerfacecolor='g', alpha = 1.0, markeredgewidth=0.0)

	upper_limit3 = np.empty(len(gmin_arr3))
	upper_limit3.fill(np.power(10.0, g_max))

	plt.fill_between(mn_arr3, gmin_arr3, upper_limit3, alpha=0.2, edgecolor='g', facecolor='k', linewidth=0.0)

	style = dict(size=15, color='g')
	plt.text(6.0, np.power(10.0, -3.2), r'$m_\nu = 0.05 \, \mathrm{eV}$', **style)

	plt.fill_between([10.0, mn_max], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

	plt.fill_between([0.0, mp], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)
	# plt.fill_betweenx([np.power(10.0, -5), 3*np.power(10.0,-4)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

	plt.fill_betweenx([np.power(10.0, -2), np.power(10.0, g_max)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)







	style = dict(size=15, color='k')
	plt.text(6.0, np.power(10.0, -3.5), r'$K^+$ decay constraint', **style)

	axes = plt.axis()
	plt.axis([0.0, mn_max, np.power(10.0, axis_min), np.power(10.0, axis_max)])

	plt.xlabel(r'$m_N / \mathrm{MeV}$')
	plt.ylabel(r'$g_\mu$')
	#plt.savefig('/Users/james/allMyStuff/Neutrinos/Constraints/plots/constraints[{},{}].pdf'.format(mp, nu_mass))
	plt.savefig('/Users/james/allMyStuff/Neutrinos/Constraints/plots/constraints[mp{}].pdf'.format(mp, nu_mass))
	print('--- Saved constraints[mp{}].pdf ---'.format(mp, nu_mass))
