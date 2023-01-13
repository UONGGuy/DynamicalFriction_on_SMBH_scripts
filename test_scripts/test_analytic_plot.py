from load_modules import *
from config import *

#run_type = 0 # indiv script selection 0 = save only / 1 = show only / 2 = save and show
fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

R_BH_CoM, V_BH_CoM, a, accel0, accel1, accel1_5 = [], [], [], [], [], []

while True:
	try:
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
		accel0 = np.append(accel0, get_snap_attr_2(fname, snap_no, 'accel_beta_0'))
		accel1 = np.append(accel1, get_snap_attr_2(fname, snap_no, 'accel_beta_1'))
		accel1_5 = np.append(accel1_5, get_snap_attr_2(fname, snap_no, 'accel_beta_1_5'))
		if snap_no == 0:
			R_BH_CoM = R_BH
			V_BH_CoM = V_BH
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
			V_BH_CoM = np.vstack((V_BH_CoM, V_BH))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1
print(accel0, accel1, accel1_5, sep='\n')

fig1 = plt.figure()
dt = np.diff(a)
dr0 = accel0[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
dr1 = accel1[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
dr1_5 = accel1_5[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
R_analyt0 = np.subtract(R_BH_CoM[:-1,3], dr0)
R_analyt0 = np.append(R_BH_CoM[0,3], R_analyt0)
R_analyt1 = np.subtract(R_BH_CoM[:-1,3], dr1)
R_analyt1 = np.append(R_BH_CoM[0,3], R_analyt1)
R_analyt1_5 = np.subtract(R_BH_CoM[:-1,3], dr1_5)
R_analyt1_5 = np.append(R_BH_CoM[0,3], R_analyt1_5)
dR = np.subtract(R_analyt1, R_analyt0)
R_analyt_corr = np.subtract(R_BH_CoM[:,3], dR)

#print(dr0[0:10], dr1[0:10], dr1_5[0:10], sep='\n')

y_max = 7.2 

plt.plot(a, R_BH_CoM[:,3], label='Numerical')
plt.plot(a, R_analyt0, label=r'Analytic $(\beta=0)$')
plt.plot(a, R_analyt1, label=r'Analytic $(\beta=1)$')
plt.plot(a, R_analyt1_5, label=r'Analytic $(\beta=1.5)$')
plt.plot(a, R_analyt_corr, label=r'Analytic corrected')
plt.vlines(decay_est_isotherm(np.linalg.norm(R_init)), ymin=0, ymax=y_max, ls='dashed', label='Isothermal decay estimate', colors='g')
plt.hlines(soft_length, xmin=0, xmax=10, ls='dashed', label='DM soft length', colors='m')
plt.hlines(soft_length_BH, xmin=0, xmax=10, ls='dashed', label='BH soft length', colors='k')

plt.xlabel(r'$t\; /(Gyr)$')
plt.xlim(0,5)
plt.ylabel(r'$R_{BH}\; /(kpc)$')
plt.ylim(0, y_max)
plt.legend(loc='upper right')

fig1.savefig('./BH_analytic_decay_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')

plt.show()
