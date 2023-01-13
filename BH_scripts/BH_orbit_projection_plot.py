from load_modules import *

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show
fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

R_BH_CoM = []

while True:
	try:
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		if snap_no == 0:
			R_BH_CoM = R_BH
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1


fig1 = plt.figure()

plt.plot(R_BH_CoM[:,0], R_BH_CoM[:,1], c='r')
plt.plot(R_BH_CoM[0,0], R_BH_CoM[0,1], 'xk', label='t = 0')
plt.plot(R_BH_CoM[-1,0], R_BH_CoM[-1,1], 'ok', label='t = 5Gyr')

#plt.title('BH orbit evolution in the X-Y plane over time \n(DM Part.=%.0e)' %(N_DM))
plt.xlabel(r'$x\; [kpc]$', fontsize=16)
plt.ylabel(r'$y\; [kpc]$', fontsize=16)

fig2 = plt.figure()

plt.plot(R_BH_CoM[:,0], R_BH_CoM[:,1], c='r')
plt.plot(R_BH_CoM[0,0], R_BH_CoM[0,1], 'xk', label='t = 0')
plt.plot(R_BH_CoM[-1,0], R_BH_CoM[-1,1], 'ok', label='t = 5Gyr')

for k in range(3):
	q = (k + 1) * 26
	plt.plot(R_BH_CoM[q,0], R_BH_CoM[q,1], 'xb')

plt.xlabel(r'$x\; [kpc]$', fontsize=16)
plt.ylabel(r'$y\; [kpc]$', fontsize=16)

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./BH_plots/BH_orbit_projection_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')
	fig2.savefig('./BH_plots/BH_wake_snap_location_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()





