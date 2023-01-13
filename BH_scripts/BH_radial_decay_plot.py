from load_modules import *
from config import *

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show
fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

R_BH_CoM = []
a = []

while True:
	try:
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
		if snap_no == 0:
			R_BH_CoM = R_BH
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1

fig1 = plt.figure()

y_max = 7.5 

plt.plot(a, R_BH_CoM[:,3], label='Numerical')
plt.vlines(decay_est_isotherm(np.linalg.norm(R_init)), ymin=0, ymax=y_max, ls='dashed', label='Isothermal decay estimate', colors='g')
plt.hlines(2.8 * soft_length, xmin=0, xmax=10, ls='dashed', label='DM soft length', colors='m')
plt.hlines(2.8 * soft_length_BH, xmin=0, xmax=10, ls='dashed', label='BH soft length', colors='k')

plt.xlabel(r'$t\; [Gyr]$', fontsize=16)
plt.xlim(0,5)
plt.ylabel(r'$R_{BH}\; [kpc]$', fontsize=16)
plt.ylim(0, y_max)
plt.legend(loc='upper right', fontsize=10)

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./BH_plots/BH_radial_decay_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()
