from load_modules import *
from config import *

#runtime = 30mins

#fig1 = plt.figure()

t_start = time.time()

snap_tot = 1
V = []
a = []
V_mean1 = np.zeros(snap_tot)
V_mean2 = np.zeros(snap_tot)
V_mean3 = np.zeros(snap_tot)
sigma2 = np.zeros(snap_tot)
local_part = 200
for k in range(snap_tot): #BH location relative to halo CoM + BH velocity for all snaps; local density + local mean particle velocity
	fname = get_snap_filename('../Output', 0)
	snap_no = 0#f'{i:03d}'
	R_DM = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
	V_DM = np.asarray(get_snap_data(fname, 1, 'Velocities'))
	R_BH = np.asarray(get_snap_data(fname, 5, 'Coordinates'))
	V_BH = np.asarray(get_snap_data(fname, 5, 'Velocities'))
	n = get_attribute(fname, 'NumPart_ThisFile')[1]
	M = np.full(n, M200 / N_DM)
	M_local = M[0:local_part]
	R_CoM = find_CoM(R_DM, M)
	a = np.append(a, get_attribute(fname, 'Time'))
	D = dist_CoM(R_CoM, R_DM)
	R_BH_CoM = np.asarray(R_BH - R_CoM)
	R_BH_CoM_mag = np.linalg.norm(R_BH_CoM)
	V_BH_mag = np.linalg.norm(V_BH)
	V_DM_mag = np.linalg.norm(V_DM, axis=1)
	V_DM_mean_3D, edges, binss = stats.binned_statistic(D[:,3], V_DM_mag, 'mean', bins=[np.linalg.norm(R_init)-0.5, np.linalg.norm(R_init)+0.5])
	V_mean1[k] = V_DM_mean_3D[0] 
	V_DM_shell_ii = np.where((D[:,3] > np.linalg.norm(R_init)-0.5) & (D[:,3] < np.linalg.norm(R_init)+0.5))
	V_mean2[k] = np.mean(np.linalg.norm(V_DM[V_DM_shell_ii], axis=1))
	if N_DM >= 1e7:
		tree = cKDTree(D[:,:-1])
		dd, ii = tree.query(R_BH_CoM, k=local_part)
		V_DM_rel = V_DM[ii]# - V_BH 
		#V_mean2[k] = np.linalg.norm(np.mean(V_DM_rel[0], axis=0))
		V_mean3[k] = np.mean(np.linalg.norm(V_DM_rel[0], axis=1))
		sigma2[k] = np.mean(np.std(V_DM_rel[0], axis=0))
#	if k==0:
#		R = np.c_[R_BH_CoM, R_BH_CoM_mag]
#		V = np.c_[V_BH, V_BH_mag]
#	else:
#		R = np.vstack((R, np.c_[R_BH_CoM, R_BH_CoM_mag]))
#		V = np.vstack((V, np.c_[V_BH, V_BH_mag]))
	V =  np.append(V, V_BH_mag)
	print(k)
	
print('snap_no =' , 82)
print('Gyr =', a)
print('Binned |v|_mean =', V_mean1)
print('|v|_mean (np.where) =', V_mean2)
print('Local_|v|_mean =', V_mean3)
print('Local_sigma =', sigma2)
print('Local_sigma/sqrt(3) =', sigma2/np.sqrt(3))
print('V_BH_mag =', V)
#with open('BH_analytic_decay_values.csv', 'w', newline='') as csvfile:
#	writer = csv.writer(csvfile)
#	writer.writerow(['a', 'R_BH_CoM', 'V_BH', 'rho_local', 'accel(0)', 'accel(1)', 'accel(1.5)'])
#	writer.writerow(a)
#	writer.writerow(R.tolist())
#	writer.writerow(V.tolist())
#	writer.writerow(rho_local)
#	writer.writerow(accel0)
#	writer.writerow(accel1)
#	writer.writerow(accel1_5)

t_finish = time.time()

print('Run time =', t_finish-t_start)

