from load_modules import *
import time

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### CODE: EXTRACT VELOCITY FOR ANG MOM ####

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

L_BH_CoM = np.cross(R_BH_CoM[:,:-1], V_BH_CoM[:,:-1], axis=1)
R_unit = R_BH_CoM[:,0:3] / R_BH_CoM[:,3][:, None]
L_mag = np.linalg.norm(L_BH_CoM, axis=1)
L_r = np.sum(L_BH_CoM * R_unit, axis=1)
L_az = np.sqrt(L_mag**2 - L_r**2)

fig_L = plt.figure()
plt.plot(a, L_mag, label='L_mag')
plt.ylabel(r'$|L_{BH}|\; [kpc\, km/s]$', fontsize=16)
plt.xlabel(r'$t\; [Gyr]$', fontsize=16)
#plt.legend(loc='upper left')

fig_L_r = plt.figure()
plt.plot(a, L_r, label='L_radial')
plt.ylabel(r'$L_{BH,r}\; [kpc\, km/s]$', fontsize=16)
plt.xlabel(r'$t\; [Gyr]$', fontsize=16)

fig_L_az = plt.figure()
plt.plot(a, L_az, label=r'$L_{BH,\phi}$')
plt.ylabel(r'$|L_{BH,\phi}|\; [kpc\, km/s]$', fontsize=16)
plt.xlabel(r'$t\; [Gyr]$', fontsize=16)

R_lim = np.argwhere(R_BH_CoM[:,3]<num_lim)
snap_no_sink = R_lim[0][0]
print('snap_no_sink =', snap_no_sink)

if snap_no_sink < snap_no:
        t_decay = get_snap_attr_2(fname, snap_no_sink, 'Time')
        L_az_res = np.mean(L_az[snap_no_sink:])
        plt.axhline(L_az_res, ls='--', label=r'$L_{BH,\phi,res}$', color='r')
        plt.axvline(t_decay, ls='--', label=r'$t_{decay}$', color='g')
        print(t_decay)

plt.legend(loc='upper right', fontsize=16)

if run_type == 0 or run_type == 2:
	fig_L.savefig('./BH_plots/BH_L_mag_DM_' + f'{N_DM:.03}' + '.pdf', dpi=600, bbox_inches='tight')
	fig_L_r.savefig('./BH_plots/BH_L_r_DM_' + f'{N_DM:.03}' + '.pdf', dpi=600, bbox_inches='tight')
	fig_L_az.savefig('./BH_plots/BH_L_az_DM_' + f'{N_DM:.03}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()


