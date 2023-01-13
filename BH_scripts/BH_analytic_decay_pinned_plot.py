from load_modules import *
from config import *

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show
fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

R_BH_CoM, V_BH_CoM, a = [], [], []
rho_local, sigma = [], []
accel0, accel1, accel1_5, accel2, accel3 = [], [], [], [], []

while True:
	try:
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
#		accel0 = np.append(accel0, get_snap_attr_2(fname, snap_no, 'accel_beta_0'))
#		accel1 = np.append(accel1, get_snap_attr_2(fname, snap_no, 'accel_beta_1'))
#		accel1_5 = np.append(accel1_5, get_snap_attr_2(fname, snap_no, 'accel_beta_1_5'))
#		accel2 = np.append(accel2, get_snap_attr_2(fname, snap_no, 'accel_beta_2'))
#		accel3 = np.append(accel3, get_snap_attr_2(fname, snap_no, 'accel_beta_3'))
		rho_local = np.append(rho_local, get_snap_attr_2(fname, snap_no, 'rho_local_pin'))
		sigma = np.append(sigma, get_snap_attr_2(fname, snap_no, 'sigma_local_pin'))
		if snap_no == 0:
			R_BH_CoM = R_BH
			V_BH_CoM = V_BH
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
			V_BH_CoM = np.vstack((V_BH_CoM, V_BH))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1
#print(accel0, accel1, accel1_5, sep='\n')

dt = np.diff(a)
#dr0 = accel0[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
#dr1 = accel1[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
#dr1_5 = accel1_5[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
dr__0 = radial_change_est(R_BH_CoM[:-1,3], np.linalg.norm(V_BH_CoM[:-1,:], axis=1), sigma[:-1], rho_local[:-1], 0) * dt * Gyr_sec
dr__1 = radial_change_est(R_BH_CoM[:-1,3], np.linalg.norm(V_BH_CoM[:-1,:], axis=1), sigma[:-1], rho_local[:-1], 1) * dt * Gyr_sec
dr__1_5 = radial_change_est(R_BH_CoM[:-1,3], np.linalg.norm(V_BH_CoM[:-1,:], axis=1), sigma[:-1], rho_local[:-1], 1.5) * dt * Gyr_sec
#dr2 = accel2[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
#dr3 = accel3[:-1] * dt * dt * (Gyr_sec)**2 / const.parsec
#R_analyt0 = np.subtract(R_BH_CoM[:-1,3], dr0)
#R_analyt0 = np.append(R_BH_CoM[0,3], R_analyt0)
#R_analyt1 = np.subtract(R_BH_CoM[:-1,3], dr1)
#R_analyt1 = np.append(R_BH_CoM[0,3], R_analyt1)
#R_analyt1_5 = np.subtract(R_BH_CoM[:-1,3], dr1_5)
#R_analyt1_5 = np.append(R_BH_CoM[0,3], R_analyt1_5)
R_analyt__0 = np.add(R_BH_CoM[:-1,3], dr__0)
R_analyt__0 = np.append(R_BH_CoM[0,3], R_analyt__0)
R_analyt__1 = np.add(R_BH_CoM[:-1,3], dr__1)
R_analyt__1 = np.append(R_BH_CoM[0,3], R_analyt__1)
R_analyt__1_5 = np.add(R_BH_CoM[:-1,3], dr__1_5)
R_analyt__1_5 = np.append(R_BH_CoM[0,3], R_analyt__1_5)
#R_analyt2 = np.subtract(R_BH_CoM[:-1,3], dr2)
#R_analyt2 = np.append(R_BH_CoM[0,3], R_analyt2)
#R_analyt3 = np.subtract(R_BH_CoM[:-1,3], dr3)
#R_analyt3 = np.append(R_BH_CoM[0,3], R_analyt3)
#dR_1 = np.subtract(R_analyt1, R_analyt0)
#dR_1_5 = np.subtract(R_analyt1_5, R_analyt0)
dR__1 = np.subtract(R_analyt__1, R_analyt__0)
dR__1_5 = np.subtract(R_analyt__1_5, R_analyt__0)
#dR_2 = np.subtract(R_analyt2, R_analyt0)
#dR_3 = np.subtract(R_analyt3, R_analyt0)
#R_analyt_corr_1 = np.subtract(R_BH_CoM[:,3], dR_1)
#R_analyt_corr_1_5 = np.subtract(R_BH_CoM[:,3], dR_1_5)
R_analyt_corr__1 = np.subtract(R_BH_CoM[:,3], dR__1)
R_analyt_corr__1_5 = np.subtract(R_BH_CoM[:,3], dR__1_5)
#R_analyt_corr_2 = np.subtract(R_BH_CoM[:,3], dR_2)
#R_analyt_corr_3 = np.subtract(R_BH_CoM[:,3], dR_3)

#print(dr0[0:10], dr1[0:10], dr1_5[0:10], sep='\n')

#### CODE: PLOTS ####

######## CODE: PLOT BETA VALUES ########

y_max = 7.2 

fig1 = plt.figure()

plt.plot(a, R_BH_CoM[:,3], label='Numerical')
#plt.plot(a, R_analyt0, label=r'Analytic $(\beta=0)$')
#plt.plot(a, R_analyt1, label=r'Analytic $(\beta=1)$')
#plt.plot(a, R_analyt1_5, label=r'Analytic $(\beta=1.5)$')
plt.plot(a, R_analyt__0,label=r'Analytic $(\beta=0)$') 
plt.plot(a, R_analyt__1, label=r'Analytic $(\beta=1)$')
plt.plot(a, R_analyt__1_5,label=r'Analytic $(\beta=1.5)$') 
#plt.plot(a, R_analyt2, label=r'Analytic $(\beta=2)$')
#plt.plot(a, R_analyt3, label=r'Analytic $(\beta=3)$')

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

plt.plot(a, R_BH_CoM[:,3], label='Numerical')
#plt.plot(a, R_analyt_corr_1, label=r'Corrected $(\beta=1)$')
#plt.plot(a, R_analyt_corr_1_5, label=r'Corrected $(\beta=1.5)$')
plt.plot(a, R_analyt_corr__1, label=r'Analytic $(\beta=1)$')
plt.plot(a, R_analyt_corr__1_5, label=r'Analytic $(\beta=1.5)$')
#plt.plot(a, R_analyt_corr_2, label=r'Corrected $(\beta=2)$')
#plt.plot(a, R_analyt_corr_3, label=r'Corrected $(\beta=3)$')

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
	fig1.savefig('./BH_plots/BH_analytic_decay_pinned_DM_' + f'{N_DM:.0e}' + '_betas.pdf', dpi=600, bbox_inches='tight')
	fig2.savefig('./BH_plots/BH_analytic_decay_pinned_DM_' + f'{N_DM:.0e}' + '_corrected.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()
