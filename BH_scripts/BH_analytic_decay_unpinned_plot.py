from load_modules import *
from config import *

#run_type = 0 # indiv script selection 0 = save only / 1 = show only / 2 = save and show
fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

R_BH_CoM, V_BH_CoM, a = [], [], []
R_BH_0, R_BH_1, R_BH_1_5 = [], [], []
V_BH_0, V_BH_1, V_BH_1_5 = [], [], []

while True:
	try:
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
		R_BH_0_hold = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH_unpin_0'))
		R_BH_1_hold = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH_unpin_1'))
		R_BH_1_5_hold = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH_unpin_1_5'))
		V_BH_0_hold = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH_unpin_0'))
		V_BH_1_hold = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH_unpin_1'))
		V_BH_1_5_hold = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH_unpin_1_5'))
		if snap_no == 0:
			R_BH_CoM = R_BH
			V_BH_CoM = V_BH
			R_BH_0 = R_BH_0_hold
			R_BH_1 = R_BH_0_hold
			R_BH_1_5 = R_BH_1_5_hold
			V_BH_0 = V_BH_0_hold
			V_BH_1 = V_BH_1_hold
			V_BH_1_5 = V_BH_1_5_hold
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
			V_BH_CoM = np.vstack((V_BH_CoM, V_BH))
			R_BH_0 = np.vstack((R_BH_0, R_BH_0_hold))
			R_BH_1 = np.vstack((R_BH_1, R_BH_1_hold))
			R_BH_1_5 = np.vstack((R_BH_1_5, R_BH_1_5_hold))
			V_BH_0 = np.vstack((V_BH_0, V_BH_0_hold))
			V_BH_1 = np.vstack((V_BH_1, V_BH_1_hold))
			V_BH_1_5 = np.vstack((V_BH_1_5, V_BH_1_5_hold))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1

dR_1 = np.subtract(R_BH_0[:,3], R_BH_1[:,3])
dR_1_5 = np.subtract(R_BH_0[:,3], R_BH_1_5[:,3])
R_analyt_1 = np.subtract(R_BH_CoM[:,3], dR_1)
R_analyt_1_5 = np.subtract(R_BH_CoM[:,3], dR_1_5)
#print(dR_1, dR_1_5,R_BH_CoM[:,3], R_BH_0, R_BH_1, R_BH_1_5, sep='\n')

#### CODE: PLOTS ####

######## CODE: PLOT BETA VALUES ########

y_max = 7.2 

fig1 = plt.figure()

plt.plot(a[:-1], R_BH_CoM[:,3], label='Numerical')
plt.plot(a[:-1], R_BH_0[:,3], label=r'Analytic $(\beta=0)$')
plt.plot(a[:-1], R_BH_1[:,3], label=r'Analytic $(\beta=1)$')
plt.plot(a[:-1], R_BH_1_5[:,3], label=r'Analytic $(\beta=1.5)$')

plt.vlines(decay_est_isotherm(np.linalg.norm(R_init)), ymin=0, ymax=y_max, ls='dashed', label='Isothermal decay estimate', colors='r')
plt.hlines(soft_length, xmin=0, xmax=10, ls='dashed', label='DM soft length', colors='m')
plt.hlines(soft_length_BH, xmin=0, xmax=10, ls='dashed', label='BH soft length', colors='k')

plt.xlabel(r'$t\; [Gyr]$')
plt.xlim(0,5)
plt.ylabel(r'$R_{BH}\; [kpc]$')
plt.ylim(0, y_max)
plt.legend(loc='upper right')

######## CODE: PLOT CORRECTED SOLUTION ########

fig2 = plt.figure()

plt.plot(a[:-1], R_BH_CoM[:,3], label='Numerical')
plt.plot(a[:-1], R_analyt_1, label=r'Analytic $(\beta=1)$')
plt.plot(a[:-1], R_analyt_1_5, label=r'Analytic $(\beta=1.5)$')

plt.vlines(decay_est_isotherm(np.linalg.norm(R_init)), ymin=0, ymax=y_max, ls='dashed', label='Isothermal decay estimate', colors='r')
plt.hlines(soft_length, xmin=0, xmax=10, ls='dashed', label='DM softening length', colors='m')
plt.hlines(soft_length_BH, xmin=0, xmax=10, ls='dashed', label='BH softening length', colors='k')

plt.xlabel(r'$t\; [Gyr]$')
plt.xlim(0,3)
plt.ylabel(r'$R_{BH}\; [kpc]$')
plt.ylim(0, y_max)
plt.legend(loc='upper right')

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./BH_plots/BH_analytic_decay_unpinned_DM_' + f'{N_DM:.0e}' + '_betas.pdf', dpi=600, bbox_inches='tight')
	fig2.savefig('./BH_plots/BH_analytic_decay_unpinned_DM_' + f'{N_DM:.0e}' + '_corrected.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()
