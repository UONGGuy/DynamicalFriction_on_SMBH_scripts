from load_modules import *
from value_paths import *

R4, R5, R6, R7, R8 = [], [], [], [], []
a4, a5, a6, a7, a8 = [], [], [], [], []
snap_tot = 180

for snap_no in range(snap_tot):

	R4_BH = np.asarray(get_snap_data_2(DM_1e4.BH_path, snap_no, 'Coordinates_BH'))
	R5_BH = np.asarray(get_snap_data_2(DM_1e5.BH_path, snap_no, 'Coordinates_BH'))
	R6_BH = np.asarray(get_snap_data_2(DM_1e6.BH_path, snap_no, 'Coordinates_BH'))
	R7_BH = np.asarray(get_snap_data_2(DM_1e7.BH_path, snap_no, 'Coordinates_BH'))
	R8_BH = np.asarray(get_snap_data_2(DM_1e8.BH_path, snap_no, 'Coordinates_BH'))

	a4 = np.append(a4, get_snap_attr_2(DM_1e4.BH_path, snap_no, 'Time'))
	a5 = np.append(a5, get_snap_attr_2(DM_1e5.BH_path, snap_no, 'Time'))
	a6 = np.append(a6, get_snap_attr_2(DM_1e6.BH_path, snap_no, 'Time'))
	a7 = np.append(a7, get_snap_attr_2(DM_1e7.BH_path, snap_no, 'Time'))
	a8 = np.append(a8, get_snap_attr_2(DM_1e8.BH_path, snap_no, 'Time'))
	
	if snap_no == 0:
		R4 = R4_BH
		R5 = R6_BH
		R6 = R6_BH
		R7 = R7_BH
		R8 = R8_BH
	else:
		R4 = np.vstack((R4, R4_BH))
		R5 = np.vstack((R5, R5_BH))
		R6 = np.vstack((R6, R6_BH))
		R7 = np.vstack((R7, R7_BH))
		R8 = np.vstack((R8, R8_BH))
			
fig1 = plt.figure()

C = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

plt.plot(a4, R4[:,3], color=C[0], label=r'$N_{DM}=10^4$')
plt.plot(a5, R5[:,3], color=C[1], label=r'$N_{DM}=10^5$')
plt.plot(a6, R6[:,3], color=C[2], label=r'$N_{DM}=10^6$')
plt.plot(a7, R7[:,3], color=C[3], label=r'$N_{DM}=10^7$')
plt.plot(a8, R8[:,3], color=C[4], label=r'$N_{DM}=10^4$')

y_max = 7.7 

plt.hlines(2.8 * 2 * DM_1e4.soft_length, xmin=0, xmax=10, ls='dashed', label=r'$2.8\,\epsilon_{BH,1e4)}$', colors=C[0])
plt.hlines(2.8 * 2 * DM_1e5.soft_length, xmin=0, xmax=10, ls='dashed', label=r'$2.8\,\epsilon_{BH,1e5)}$', colors=C[1])
plt.hlines(2.8 * 2 * DM_1e6.soft_length, xmin=0, xmax=10, ls='dashed', label=r'$2.8\,\epsilon_{BH,1e6)}$', colors=C[2])
plt.hlines(2.8 * 2 * DM_1e7.soft_length, xmin=0, xmax=10, ls='dashed', label=r'$2.8\,\epsilon_{BH,1e7)}$', colors=C[3])
plt.hlines(2.8 * 2 * DM_1e8.soft_length, xmin=0, xmax=10, ls='dashed', label=r'$2.8\,\epsilon_{BH,1e8)}$', colors=C[4])

plt.vlines(decay_est_isotherm(np.linalg.norm(R_init)), ymin=0, ymax=y_max, ls='dashed', label=r'$t_{decay,isotherm}$', colors=C[5])

plt.xlabel(r'$t\; [Gyr]$', fontsize=14)
plt.xlim(0,4.5)
plt.ylabel(r'$R_{BH}\; [kpc]$', fontsize=14)
plt.ylim(0, y_max)
plt.legend(loc='upper right', fontsize=12)

fig1.savefig('./BH_plots_compare/BH_decay_comparison.pdf', dpi=600, bbox_inches='tight')

plt.grid()

plt.show()
