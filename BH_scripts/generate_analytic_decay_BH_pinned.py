from load_modules import *
from config import *

#runtime() use sintr -A DIRAC-DP012-SL4-CPU -N16 -n896 -t 2:00:00 -p cclake#

fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

ti = time.time()

R_BH_CoM = []
V_BH_CoM = []
a = []
rho_local = []
sigma = []
n_part_local = 64
local_radius = []
local_part_rel = []

while True:
	try:
		R_DM = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_DM'))
		V_DM = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_DM'))
		R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
		V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
		n = get_snap_attr_2(fname, snap_no, 'NumPart_DM')
		M = np.full(n, M200 / N_DM)
		M_local = M[0:n_part_local]
		tree = cKDTree(R_DM[:,:-1])
		dd, ii = tree.query(R_BH[0,:-1], k=n_part_local)
		V_DM_local = V_DM[ii] 
		v_thresh_ii = np.where(np.linalg.norm(V_DM_local, axis=1) < np.linalg.norm(V_BH))
		local_part_rel = np.append(local_part_rel, len(v_thresh_ii[0]))
		local_radius =  np.append(local_radius, max(dd))
#		V_DM_local_thresh = V_DM_local[v_thresh_ii]
#		sigma_3 = np.std(V_DM_local, axis=0)#_thresh, axis=0)
		sigma = np.append(sigma, np.std(V_DM_local))#np.average(sigma_3))
		if snap_no == 0:
			R_BH_CoM = R_BH
			V_BH_CoM = np.c_[V_BH, np.linalg.norm(V_BH)]
		else:
			R_BH_CoM = np.vstack((R_BH_CoM, R_BH))
			V_BH_CoM = np.vstack((V_BH_CoM, np.c_[V_BH, np.linalg.norm(V_BH)]))
		print(snap_no)
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		tf = time.time()
		print(f'Error occured during snap_{snap_no:3d}')
		print('Runtime = ', tf-ti)
		break
	else:
		if np.isnan(sigma[snap_no]):
			print('snap_no = ', snap_no)
			print('V_BH', V_BH, 'tree', V_DM_local, 'local_part_rel', local_part_rel[snap_no], 'sigma_3,', sigma_3, 'sigma', sigma, sep='\n')	
		snap_no += 1

rho_local = [local_part_rel[k] * M200 / N_DM / (4 / 3 * np.pi * local_radius[k]**3) for k in range(len(a))]
#accel0 = [accel_est_A(R_BH_CoM[k,3], V_BH_CoM[k,3], sigma[k], rho_local[k], 0) for k in range(len(a))]
#accel1 = [accel_est_A(R_BH_CoM[k,3], V_BH_CoM[k,3], sigma[k], rho_local[k], 1) for k in range(len(a))]
#accel1_5 = [accel_est_A(R_BH_CoM[k,3], V_BH_CoM[k,3], sigma[k], rho_local[k], 1.5) for k in range(len(a))]
#accel2 = [accel_est_A(R_BH_CoM[k,3], V_BH_CoM[k,3], sigma[k], rho_local[k], 2) for k in range(len(a))]
#accel3 = [accel_est_A(R_BH_CoM[k,3], V_BH_CoM[k,3], sigma[k], rho_local[k], 3) for k in range(len(a))]

#print('sigma', sigma, 'local_radius', local_radius, 'local_part_rel', local_part_rel, 'rho_local', rho_local, sep='\n')

hf = h5py.File('values_BH_DM_' + f'{N_DM:.0e}' + '.h5', 'r+')
for i in range(len(a)):
	snap_i = hf[f'snap_{i:03d}']
	snap_i.attrs['rho_local_pin'] = rho_local[i]
	snap_i.attrs['sigma_local_pin'] = sigma[i]
#	snap_i.attrs['accel_beta_0'] = accel0[i]
#	snap_i.attrs['accel_beta_1'] = accel1[i]
#	snap_i.attrs['accel_beta_1_5'] = accel1_5[i]
#	snap_i.attrs['accel_beta_2'] = accel2[i]
#	snap_i.attrs['accel_beta_3'] = accel3[i]

tf = time.time()

print('Runtime =', tf-ti)


