import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ja')

txs_df = pd.read_csv(filepath_or_buffer='txsEvents.csv',
                     header=0, sep='   ', engine='python')

nontxs_df = pd.read_csv(filepath_or_buffer='nontxsEvents.csv',
                        header=0, sep='   ', engine='python')

mjd_arr = np.array(txs_df['MJD'])
nontxs_mjd_arr = np.array(nontxs_df['MJD'])
logE_arr = np.array(txs_df['log10(Ereco)'])
nontxs_logE_arr = np.array(nontxs_df['log10(Ereco)'])

Emu_TeV = 10**(-3) * 10**logE_arr
Enu_TeV = 1.92 * Emu_TeV**(1.14)

# # Data Bar Chart to show individual Events
# plt.figure()
# plt.bar(nontxs_mjd_arr, nontxs_logE_arr,
#         facecolor='k',
#         width=1.0,
#         linewidth=0,
#         align='edge')
# plt.bar(mjd_arr, logE_arr,
#         facecolor='r',
#         width=1.0,
#         linewidth=0,
#         align='edge')
# # plt.plot(mjd_arr, logE_arr, 'r-.', linewidth=0.5)
# plt.xlabel('Modified Julian Date')
# plt.ylabel(r'$\log_{10} \left[ \frac{\hat{E}_{\mu}}{\textrm{\tiny GeV}} \right]$')
# plt.savefig('discover1.pdf')
#
# # Histogram plot of Reconstructed Neutrino Energies
# plt.figure()
# plt.subplot(2, 1, 1)
# plt.hist(Emu_TeV, 100, facecolor='g', alpha=0.1, cumulative=True, linewidth=0, label='Cumulative')
# plt.hist(Emu_TeV, 20, facecolor='r', alpha=0.5)
# plt.legend(loc='best')
# axes = plt.axis()
# plt.xlabel(r'$\hat{E}_{\mu} / \textrm{TeV}$')
# plt.ylabel(r'Number of Events')
# new_axes = [0, max(Emu_TeV), axes[2], axes[3]]
# plt.axis(new_axes)
# plt.subplot(2, 1, 2)
# plt.hist(Enu_TeV, 100, facecolor='g', alpha=0.1, cumulative=True, linewidth=0)
# plt.hist(Enu_TeV, 20, facecolor='b', alpha=0.5, cumulative=False)
# plt.xlabel(r'$E_{\nu} / \textrm{TeV}$')
# plt.ylabel(r'Number of Events')
# axes = plt.axis()
# new_axes = [0, max(Enu_TeV), axes[2], axes[3]]
# plt.axis(new_axes)
# plt.savefig('discover2.pdf')
