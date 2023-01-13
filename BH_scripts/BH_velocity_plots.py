from load_modules import *
import time

#run_type = 1 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### CODE: EXTRACT VELOCITY ####

fname = f'values_BH_DM_{N_DM:.0e}.h5'
num_lim = 2.8 * soft_length_BH
snap_no = 0
R_BH_CoM = []
V_BH_CoM = []
a = []

while True:
	try:
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
		if snap_no == 0:
			R_BH_CoM = R_BH
			V_BH_CoM = V_BH
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
			V_BH_CoM = np.vstack((V_BH_CoM, V_BH))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1

######## CODE: ROC PLOTS ########

R_unit = R_BH_CoM[:,0:3] / R_BH_CoM[:,3][:, None]
V_mag = np.linalg.norm(V_BH_CoM, axis=1)
V_r = np.sum(V_BH_CoM * R_unit, axis=1)
V_az = np.sqrt(V_mag**2 - V_r**2)

fig_V = plt.figure()
plt.plot(a, V_mag, label='V_mag')
plt.ylabel(r'$|v_{BH}|\; [km/s]$')
plt.xlabel(r'$t\; [Gyr]$')
#plt.legend(loc='upper left')

fig_V_r = plt.figure()
plt.plot(a, V_r, label='V_radial')
plt.ylabel(r'$v_{BH,r}\; [km/s]$')
plt.xlabel(r'$t\; [Gyr]$')

fig_V_az = plt.figure()
plt.plot(a, V_az, label=r'$|V_{BH,\phi}|$')
plt.ylabel(r'$|v_{BH,\phi}|\; [km/s]$', fontsize=18)
plt.xlabel(r'$t\; [Gyr]$', fontsize=18)

######## CODE: RESIDUAL VELOCITY ########

R_lim = np.argwhere(R_BH_CoM[:,3]<num_lim)
snap_no_sink = R_lim[0][0]
print('snap_no_sink =', snap_no_sink)

if snap_no_sink < snap_no:
	t_decay = get_snap_attr_2(fname, snap_no_sink, 'Time')
	V_az_res = np.mean(V_az[snap_no_sink:])
	plt.axhline(V_az_res, ls='--', label=r'$V_{BH,\phi,res}$', color='r')
	plt.axvline(t_decay, ls='--', label=r'$t_{decay}$', color='g')
#	plt.hlines(V_az_res, xmin=0, xmax=5, ls='--', label=r'$V_{BH,\phi,res}$', color='r')
#	plt.vlines(t_decay, ymin=35, ymax=56, ls='--', label=r'$t_{decay}$', color='g')
	print(t_decay)

plt.legend(loc='upper right', fontsize=16)


if run_type == 0 or run_type == 2:
	fig_V.savefig('./BH_plots/BH_V_mag_DM_' + f'{N_DM:.03}' + '.pdf', dpi=600, bbox_inches='tight')
	fig_V_r.savefig('./BH_plots/BH_V_r_DM_' + f'{N_DM:.03}' + '.pdf', dpi=600, bbox_inches='tight')
	fig_V_az.savefig('./BH_plots/BH_V_az_DM_' + f'{N_DM:.03}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()


