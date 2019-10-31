import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import sys
warnings.filterwarnings("ignore")
#plt.style.use('ja')
data_dir = 'final_data/'

if __name__ == '__main__':
	mp = float(sys.argv[1]) # MeV
	min_nu_mass = float(sys.argv[2]) # eV

	hierarchy = 'normal'
	filename = 'mp' + str(mp) + 'numin' + str(min_nu_mass) + hierarchy.upper() + '.csv'
	data = pd.read_csv(data_dir + filename, index_col=0)
	print('--- Opened ' + data_dir + filename + ' ---')
	gmin_arrN = data['gmin_arr']
	mn_arrN = data['mn_arr']


	hierarchy = 'inverted'
	filename = 'mp' + str(mp) + 'numin' + str(min_nu_mass) + hierarchy.upper() + '.csv'
	data = pd.read_csv(data_dir + filename, index_col=0)
	print('--- Opened ' + data_dir + filename + ' ---')
	gmin_arrI = data['gmin_arr']
	mn_arrI = data['mn_arr']

	mn_max = 10.0 # MeV
	g_min = -5.0
	g_max = -1.0
	axis_min = -4.0
	axis_max = -1.0

	fs = 16

	plt.figure()
	#plt.rc('axes', linewidth=2.0)

	plt.plot([mp, mp], [np.power(10.0, -4), np.power(10.0, -2)], c='k', linewidth=1.5)
	plt.plot([mp, 10.0], [np.power(10.0, -2), np.power(10.0, -2)], c='k', linewidth=1.5)
	#plt.plot([10.0, 10.0], [np.power(10.0, -2), np.power(10.0, -4)], c='k', linewidth=1.0)

	style = dict(size=fs, color='k')
	#plt.text(0.5, np.power(10.0, -1.5), r'$m_\delta = {} \,$'.format(mp) + r'$\mathrm{MeV}$, ' + r'$m^{\mathrm{min}}_\nu$' + r'$= {} \,$'.format(min_nu_mass) + r'$\mathrm{eV}$', **style, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
	plt.text(6.0, np.power(10.0, -3.7), r'$m_\delta = {} \,$'.format(mp) + r'$\mathrm{MeV}$', **style)
	plt.text(6.0, np.power(10.0, -3.9), r'$m^{\mathrm{min}}_\nu$' + r'$= {} \,$'.format(min_nu_mass) + r'$\mathrm{eV}$', **style)

	normal_color = '#357DED'
	inv_color = '#D1495B'

	plt.semilogy(mn_arrN, gmin_arrN, linestyle='-', color=normal_color, label='Normal Hierarchy')

	upper_limit1 = np.empty(len(gmin_arrN))
	upper_limit1.fill(np.power(10.0, g_max))

	plt.fill_between(mn_arrN, gmin_arrN, upper_limit1, alpha=0.1, edgecolor='k', facecolor=normal_color, linewidth=2.0)

	style = dict(size=fs, color=normal_color)
	#plt.text(2.6, np.power(10.0, -3.5), 'Normal Hierarchy', **style)

	plt.semilogy(mn_arrI, gmin_arrI, linestyle='-', color=inv_color, label='Inverted Hierarchy')

	upper_limit2 = np.empty(len(gmin_arrI))
	upper_limit2.fill(np.power(10.0, g_max))

	plt.fill_between(mn_arrI, gmin_arrI, upper_limit2, alpha=0.1, edgecolor='k', facecolor=inv_color, linewidth=0.0)

	style = dict(size=fs, color=inv_color)
	#plt.text(2.6, np.power(10.0, -3.7), 'Inverted Hierarchy', **style)

	plt.fill_between([10.0, mn_max], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)

	plt.fill_between([0.0, mp], [np.power(10.0, g_min), np.power(10.0, g_min)], [np.power(10.0, g_max), np.power(10.0, g_max)], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)


	plt.fill_betweenx([np.power(10.0, -2), np.power(10.0, g_max)], [mp, mp], [10.0, 10.0], alpha=0.1, edgecolor='k', facecolor='k', linewidth=0.0)




	style = dict(size=fs, color='k')
	plt.text(mp + 0.1, np.power(10.0, -1.95), r'$K^+$ decay', fontdict=style)
	plt.text(mp + 0.1, np.power(10.0, -3.4), r'$m_\delta < m_N$', fontdict=style, rotation=90.0)
	# plt.annotate('', xy=(10.0, np.power(10.0, -3.8)), xytext=(8.0, np.power(10.0, -3.55)),
	#             arrowprops=dict(arrowstyle="->",
	#                             connectionstyle="arc3,rad=0.2"),
	#             )

	axes = plt.axis()
	plt.axis([0.0, mn_max, np.power(10.0, axis_min), np.power(10.0, axis_max)])
	plt.xticks(np.arange(0, 11, 1))
	plt.gca().set_xticks(np.arange(0, 10, 0.2), minor=True)
	ax = plt.gca()
	ax.tick_params(which='minor', length=1.5)

	plt.xlabel(r'$m_N \, [\mathrm{MeV}]$')
	plt.ylabel(r'$g_\mu$')
	plt.legend(fontsize=12, loc='upper left')

	plt.savefig('final_plots/constraints[mp{}numin{}].pdf'.format(mp, min_nu_mass))
	print('--- Saved constraints[mp{}numin{}].pdf ---'.format(mp, min_nu_mass))
