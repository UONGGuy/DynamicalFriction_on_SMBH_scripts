from load_modules import *

run_type = 1 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### CODE: EXTRACT VALUES ####

fig1 = plt.figure()

R = []
a = []
local_part = 64
for k in range(201): #BH location relative to halo CoM + BH velocity for all snaps; local density + local mean particle velocity
	fname = get_snap_filename('../Output', k)
	snap_no = f'{i:03d}'
	R_DM = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
	V_DM = np.asarray(get_snap_data(fname, 1, 'Velocities'))
	R_BH = np.asarray(get_snap_data(fname, 5, 'Coordinates'))
	V_BH = np.asarray(get_snap_data(fname, 5, 'Velocities'))
	n = get_attribute(fname, 'NumPart_ThisFile')[1]
	M = np.full(n,M200 /N_DM)
	M_local = M[0:local_part]
	R_CoM = find_CoM(R_DM, M)
	a = np.append(a, get_attribute(fname, 'Time'))
	D = dist_CoM(R_CoM, R_DM)
	R_BH_CoM = np.asarray(R_BH - R_CoM)
	R_BH_CoM_mag = np.linalg.norm(R_BH_CoM)
	V_BH_CoM_mag = np.linalg.norm(V_BH)
#	if N_DM >= 1e7:
#		tree = cKDTree(D[:,:-1])
#		dd, ii = tree.query(R_BH_CoM, k=local_part)
##		V_local_CoM = find_CoM([V_DM[index] for index in ii][0], M_local)
##		V_DM_local_rel = [V_DM[index] - V_local_CoM for index in ii][0]
#		V_DM_rel = [V_DM[index] - V_BH for index in ii][0]
#		v_thresh_ii = np.where(np.linalg.norm(v_DM_rel, axis=1) < np.linalg.norm(V_BH))
#		rho_local = len(v_thresh_ii) * M200 / N_DM / (4 / 3 * np.pi * max(dd[0])
	if k==0:
		R = np.c_[R_BH_CoM, R_BH_CoM_mag]
		V = np.c_[V_BH, V_BH_CoM_mag]
#		print(dd[0][-1])
	else:
		R = np.vstack((R, np.c_[R_BH_CoM, R_BH_CoM_mag]))
		V = np.vstack((V, np.c_[V_BH, V_BH_CoM_mag]))

with open('IC_BH_velocity.csv', 'r') as csvfile:
	csvreader = csv.reader(csvfile)
	header = next(csvreader)
	values = next(csvreader)

values = [np.float32(j) for j in values]

######## CODE: RADIAL ORBIT + DECAY EST ########

t_decay = values[header.index('t_decay [Gyr]')]
t_decay_isotherm = values[header.index('t_decay_isothermal [Gyr]')]

y_max = 9.5

plt.plot(a, R[:,3], label='BH distance from halo centre', c='r')
plt.vlines(t_decay, ymin=0, ymax=y_max, ls='dashed', label='Analytic decay estimate', colors='b')
plt.vlines(t_decay_isotherm, ymin=0, ymax=y_max, ls='dashed', label='Isothermal decay estimate', colors='g')
plt.hlines(soft_length, xmin=0, xmax=10, ls='dashed', label='Halo soft length', colors='m')
plt.hlines(soft_length_BH, xmin=0, xmax=10, ls='dashed', label='BH soft length', colors='k')

plt.title('BH orbit evolution in time \n(DM Part.=%.0e)' %(N_DM))
plt.xlabel('Time [Gyr]')
plt.xlim(0, 5)
plt.ylabel('Distance from CoM [kpc]')
plt.ylim(0, y_max) #adjust if necessary/make fn involving ceil to relate to r_init
plt.legend(loc='upper right')

######## CODE: ORBIT TRAJECTORY PROJECTION (XY) ########

fig2 = plt.figure()

plt.plot(R[:,0], R[:,1], label='BH distance from halo centre', c='r')
plt.plot(R[0,0], R[0,1], 'xk', label='t = 0')
plt.plot(R[-1,0], R[-1,1], 'ok', label='t = 5Gyr')

plt.title('Projection of BH orbit in the X-Y plane over time \n(DM Part.=%.0e)' %(N_DM))
plt.xlabel('x-distance from CoM [kpc]')
plt.ylabel('y-distance from CoM [kpc]')
plt.legend(loc='upper right')

######## CODE: VELOCITY ROC ########

fig3 = plt.figure()
dV = np.diff(V[:,3])
dt = np.diff(a)
dV_dt = np.divide(dV, dt)

plt.plot(a[:-1], dV_dt, label='BH distance from halo centre', c='r')
plt.plot(a[0], dV_dt[0], 'xk', label='start')
plt.plot(a[-2], dV_dt[-1], 'ok', label='end')
plt.title('Rate of change of BH velocity \n(DM Part.=%.0e)' %(N_DM))
plt.xlabel('Time [Gyr]')
plt.ylabel('Change in velocity [km/s]')
plt.legend(loc='upper right')

######## CODE: ANGULAR MOMENTUM ROC ########

fig4 = plt.figure()
L = np.cross(R[:,0:3], V[:,0:3])
L_mag = np.linalg.norm(L, axis=1)
dL = np.diff(L_mag)
dL_dt = np.divide(dL, dt)

plt.plot(a[:-1], dL_dt, label='BH distance from halo centre', c='r')
plt.plot(a[0], dL_dt[0], 'xk', label='start')
plt.plot(a[-2], dL_dt[-1], 'ok', label='end')
plt.title('Rate of change of BH angular momentum \n(DM Part.=%.0e)' %(N_DM))
plt.xlabel('Time [Gyr]')
plt.ylabel('Change in angular momentum [kpc km/s]')
plt.legend(loc='upper right')

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./BH_plots/BH_orbit_radial_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')
	fig2.savefig('./BH_plots/BH_orbit_projection_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')
	fig3.savefig('./BH_plots/BH_RoC_Vel_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')
	fig4.savefig('./BH_plots/BH_RoC_AngMom_DM' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()
