from load_modules import *
from config import *

#keep as generator/plotter separate -> create separate file for extra values here

ti = time.time()
fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

R_BH_unpin_0 = np.asarray(get_snap_data_2(fname, 0, 'Coordinates_BH'))
V_BH_unpin_0 = np.asarray(get_snap_data_2(fname, 0, 'Velocities_BH'))
R_BH_unpin_1 = np.asarray(get_snap_data_2(fname, 0, 'Coordinates_BH'))
V_BH_unpin_1 = np.asarray(get_snap_data_2(fname, 0, 'Velocities_BH'))
R_BH_unpin_1_5 = np.asarray(get_snap_data_2(fname, 0, 'Coordinates_BH'))
V_BH_unpin_1_5 = np.asarray(get_snap_data_2(fname, 0, 'Velocities_BH'))
a = []
rho_local_unpin_0 = []
sigma_local_unpin_0 = []
rho_local_unpin_1 = []
sigma_local_unpin_1 = []
rho_local_unpin_1_5 = []
sigma_local_unpin_1_5 = []
dr_0, dr_1, dr_1_5 = [], [], []
dv_0, dv_1, dv_1_5 = [], [], []
n_part_local = 64
rot_axis = np.array([0, 0, 1])

while True:
	try:
		R_DM = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_DM'))
		V_DM = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_DM'))
		R_BH_0 = np.asarray([R_BH_unpin_0[snap_no,:]])
		V_BH_0 = np.asarray([V_BH_unpin_0[snap_no,:]])
		R_BH_1 = np.asarray([R_BH_unpin_1[snap_no,:]])
		V_BH_1 = np.asarray([V_BH_unpin_1[snap_no,:]])
		R_BH_1_5 = np.asarray([R_BH_unpin_1_5[snap_no,:]])
		V_BH_1_5 = np.asarray([V_BH_unpin_1_5[snap_no,:]])
		V_BH_0_mag = np.linalg.norm(V_BH_0)
		V_BH_1_mag = np.linalg.norm(V_BH_1)
		V_BH_1_5_mag = np.linalg.norm(V_BH_1_5)
		a = np.append(a, get_snap_attr_2(fname, snap_no, 'Time'))
		print('Snap', snap_no, 'Time =', a[snap_no], 'Gyr')
		dt = get_snap_attr_2(fname, snap_no + 1, 'Time') - a[snap_no]
		n = get_snap_attr_2(fname, snap_no, 'NumPart_DM')
		M = np.full(n, M200 / N_DM)
#LOCAL DENSITY AND DISPERSION
		M_local = M[0:n_part_local]
		tree = cKDTree(R_DM[:,:-1])
	#beta=0
		dd_0, ii_0 = tree.query(R_BH_0[0,:-1], k=n_part_local)
		V_DM_local_0 = V_DM[ii_0] 
		v_thresh_ii_0 = np.where(np.linalg.norm(V_DM_local_0, axis=1) < np.linalg.norm(V_BH_0))
		local_part_rel_0 = len(v_thresh_ii_0[0])
		local_radius_0 =  max(dd_0)
		sigma_local_unpin_0 = np.append(sigma_local_unpin_0, np.std(V_DM_local_0))
		rho_local_unpin_0 = np.append(rho_local_unpin_0, local_part_rel_0 * M200 / N_DM / (4 / 3 * np.pi * local_radius_0**3))
	#beta=1	
		dd_1, ii_1 = tree.query(R_BH_1[0,:-1], k=n_part_local)
		V_DM_local_1 = V_DM[ii_1] 
#		print('V_DM_local_1', V_DM_local_1)
		v_thresh_ii_1 = np.where(np.linalg.norm(V_DM_local_1, axis=1) < np.linalg.norm(V_BH_1))
		local_part_rel_1 = len(v_thresh_ii_1[0])
		local_radius_1 =  max(dd_1)
		sigma_local_unpin_1 = np.append(sigma_local_unpin_1, np.std(V_DM_local_1))
		rho_local_unpin_1 = np.append(rho_local_unpin_1, local_part_rel_1 * M200 / N_DM / (4 / 3 * np.pi * local_radius_1**3))
	#beta=1.5
		dd_1_5, ii_1_5 = tree.query(R_BH_1_5[0,:-1], k=n_part_local)
		V_DM_local_1_5 = V_DM[ii_1_5] 
		v_thresh_ii_1_5 = np.where(np.linalg.norm(V_DM_local_1_5, axis=1) < np.linalg.norm(V_BH_1_5))
		local_part_rel_1_5 = len(v_thresh_ii_1_5[0])
		local_radius_1_5 =  max(dd_1_5)
		sigma_local_unpin_1_5 = np.append(sigma_local_unpin_1_5, np.std(V_DM_local_1_5))
		rho_local_unpin_1_5 = np.append(rho_local_unpin_1_5, local_part_rel_1_5 * M200 / N_DM / (4 / 3 * np.pi * local_radius_1_5**3))
#CHANGE IN RADIAL DISTANCE [kpc]
		dr_0 = np.append(dr_0, radial_change_est(R_BH_0[0,3], V_BH_0_mag, sigma_local_unpin_0[snap_no], rho_local_unpin_0[snap_no], 0) * dt * Gyr_sec)
		dr_1 = np.append(dr_1, radial_change_est(R_BH_1[0,3], V_BH_1_mag, sigma_local_unpin_1[snap_no], rho_local_unpin_1[snap_no], 1) * dt * Gyr_sec)
		dr_1_5 = np.append(dr_1_5, radial_change_est(R_BH_1_5[0,3], V_BH_1_5_mag, sigma_local_unpin_1_5[snap_no], rho_local_unpin_1_5[snap_no], 1.5) * dt * Gyr_sec)
		print('dr_0', dr_0[snap_no])
		print('dr_1', dr_1[snap_no])
		print('dr_1_5', dr_1_5[snap_no])
#CHANGE IN VELOCITY [kpc s^{-1}]
		dv_0 = np.append(dv_0, velocity_change_est(R_BH_0[0,3], V_BH_0_mag, sigma_local_unpin_0[snap_no], rho_local_unpin_0[snap_no], 0) * dt * Gyr_sec)
		dv_1 = np.append(dv_1, velocity_change_est(R_BH_1[0,3], V_BH_1_mag, sigma_local_unpin_1[snap_no], rho_local_unpin_1[snap_no], 1) * dt * Gyr_sec)
		dv_1_5 = np.append(dv_1_5, velocity_change_est(R_BH_1_5[0,3], V_BH_1_5_mag, sigma_local_unpin_1_5[snap_no], rho_local_unpin_1_5[snap_no], 1.5) * dt * Gyr_sec)
#CHANGE IN ANGLE [radians]
		d_theta_0 = V_BH_0_mag / R_BH_0[0,3] * dt * Gyr_sec / const.parsec
		d_theta_1 = V_BH_1_mag / R_BH_1[0,3] * dt * Gyr_sec / const.parsec
		d_theta_1_5 = V_BH_1_5_mag / R_BH_1_5[0,3] * dt * Gyr_sec / const.parsec
#ROTATE R_BH
		rot_vec_0 = d_theta_0 * rot_axis
		rot_vec_1 = d_theta_1 * rot_axis
		rot_vec_1_5 = d_theta_1_5 * rot_axis
		rotation_0 = Rot.from_rotvec(rot_vec_0)
		rotation_1 = Rot.from_rotvec(rot_vec_1)
		rotation_1_5 = Rot.from_rotvec(rot_vec_1_5)
		R_BH_rot_0 = np.asarray([rotation_0.apply(R_BH_0[0,:-1])])
		R_BH_rot_0_mag = np.linalg.norm(R_BH_rot_0)
		R_BH_rot_1 = np.asarray([rotation_1.apply(R_BH_1[0,:-1])])
		R_BH_rot_1_mag = np.linalg.norm(R_BH_rot_1)
		R_BH_rot_1_5 = np.asarray([rotation_1_5.apply(R_BH_1_5[0,:-1])])
		R_BH_rot_1_5_mag = np.linalg.norm(R_BH_rot_1_5)
		print('R_BH_rot_0', R_BH_rot_0)
		print('R_BH_rot_1', R_BH_rot_1)
		print('R_BH_rot_1_5', R_BH_rot_1_5)
#ROTATE V_BH
		V_BH_rot_0 = np.asarray([rotation_0.apply(V_BH_0[0])])
		V_BH_rot_0_mag = np.linalg.norm(V_BH_rot_0)
		V_BH_rot_1 = np.asarray([rotation_1.apply(V_BH_1[0])])
		V_BH_rot_1_mag = np.linalg.norm(V_BH_rot_1)
		V_BH_rot_1_5 = np.asarray([rotation_1.apply(V_BH_1_5[0])])
		V_BH_rot_1_5_mag = np.linalg.norm(V_BH_rot_1_5)
#COMPUTE NEW R_BH
		R_BH_rot_0_unit = R_BH_rot_0 / R_BH_rot_0_mag
		R_BH_new_0 = R_BH_rot_0 + dr_0[snap_no] * R_BH_rot_0_unit
		R_BH_rot_1_unit = R_BH_rot_1 / R_BH_rot_1_mag
		R_BH_new_1 = R_BH_rot_1 + dr_1[snap_no] * R_BH_rot_1_unit
		R_BH_rot_1_5_unit = R_BH_rot_1_5 / R_BH_rot_1_5_mag
		R_BH_new_1_5 = R_BH_rot_1_5 + dr_1_5[snap_no] * R_BH_rot_1_5_unit
#COMPUTE NEW V_BH
		V_BH_rot_0_unit = V_BH_rot_0 / V_BH_rot_0_mag
		V_BH_new_0 = V_BH_rot_0 + dv_0[snap_no] * V_BH_rot_0_unit
		V_BH_rot_1_unit = V_BH_rot_1 / V_BH_rot_1_mag
		V_BH_new_1 = V_BH_rot_1 + dv_1[snap_no] * V_BH_rot_1_unit
		V_BH_rot_1_5_unit = V_BH_rot_1_5 / V_BH_rot_1_5_mag
		V_BH_new_1_5 = V_BH_rot_1_5 + dv_1_5[snap_no] * V_BH_rot_1_5_unit
#UPDATE R_BH_unpin, V_BH_unpin
		R_BH_unpin_0 = np.vstack((R_BH_unpin_0, np.c_[R_BH_new_0, np.linalg.norm(R_BH_new_0)]))
		R_BH_unpin_1 = np.vstack((R_BH_unpin_1, np.c_[R_BH_new_1, np.linalg.norm(R_BH_new_1)]))
		R_BH_unpin_1_5 = np.vstack((R_BH_unpin_1_5, np.c_[R_BH_new_1_5, np.linalg.norm(R_BH_new_1_5)]))
		V_BH_unpin_0 = np.vstack((V_BH_unpin_0, V_BH_new_0))
		V_BH_unpin_1 = np.vstack((V_BH_unpin_1, V_BH_new_1))
		V_BH_unpin_1_5 = np.vstack((V_BH_unpin_1_5, V_BH_new_1_5))
				
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		tf = time.time()
		print(f'Error on snap_{snap_no:3d}')
		print('Runtime = ', tf-ti)
		break
	else:
		#if np.isnan(sigma[snap_no]):
		#	print('snap_no = ', snap_no)
		#	print('V_BH', V_BH, 'tree', V_DM_local, 'local_part_rel', local_part_rel[snap_no], 'sigma_3,', sigma_3, 'sigma', sigma, sep='\n')	
		snap_no += 1

print('R_BH_unpin_0', R_BH_unpin_0, sep='\n')
print('R_BH_unpin_1', R_BH_unpin_1, sep='\n')
print('R_BH_unpin_1_5', R_BH_unpin_1_5, sep='\n')

#### CODE: SAVE VALUES TO HDF5 FILE ####

hf = h5py.File('values_BH_DM_' + f'{N_DM:.0e}' + '.h5', 'r+')

lst_keys = list(hf['snap_000'].keys())

if not 'Coordinates_BH_unpin_1' in lst_keys:
	for i in range(len(a)):
		snap_i = hf[f'snap_{i:03d}']
		snap_i.attrs['rho_local_unpin_0'] = rho_local_unpin_0[i]
		snap_i.attrs['sigma_local_unpin_0'] = sigma_local_unpin_0[i]
		snap_i.attrs['rho_local_unpin_1'] = rho_local_unpin_1[i]
		snap_i.attrs['sigma_local_unpin_1'] = sigma_local_unpin_1[i]
		snap_i.attrs['rho_local_unpin_1_5'] = rho_local_unpin_1_5[i]
		snap_i.attrs['sigma_local_unpin_1_5'] = sigma_local_unpin_1_5[i]
		snap_i.create_dataset('Coordinates_BH_unpin_0', data=np.asarray([R_BH_unpin_0[i]]))
		snap_i.create_dataset('Coordinates_BH_unpin_1', data=np.asarray([R_BH_unpin_1[i]]))
		snap_i.create_dataset('Coordinates_BH_unpin_1_5', data=np.asarray([R_BH_unpin_1_5[i]]))
		snap_i.create_dataset('Velocities_BH_unpin_0', data=np.asarray([V_BH_unpin_0[i]]))
		snap_i.create_dataset('Velocities_BH_unpin_1', data=np.asarray([V_BH_unpin_1[i]]))
		snap_i.create_dataset('Velocities_BH_unpin_1_5', data=np.asarray([V_BH_unpin_1_5[i]]))

if 'Coordinates_BH_unpin_1' in lst_keys:
	for i in range(len(a)):
		snap_i = hf[f'snap_{i:03d}']
		snap_i.attrs['rho_local_unpin_0'] = rho_local_unpin_0[i]
		snap_i.attrs['sigma_local_unpin_0'] = sigma_local_unpin_0[i]
		snap_i.attrs['rho_local_unpin_1'] = rho_local_unpin_1[i]
		snap_i.attrs['sigma_local_unpin_1'] = sigma_local_unpin_1[i]
		snap_i.attrs['rho_local_unpin_1_5'] = rho_local_unpin_1_5[i]
		snap_i.attrs['sigma_local_unpin_1_5'] = sigma_local_unpin_1_5[i]
		del snap_i['Coordinates_BH_unpin_0']
		del snap_i['Coordinates_BH_unpin_1']
		del snap_i['Coordinates_BH_unpin_1_5']
		del snap_i['Velocities_BH_unpin_0']
		del snap_i['Velocities_BH_unpin_1']
		del snap_i['Velocities_BH_unpin_1_5']
		snap_i.create_dataset('Coordinates_BH_unpin_0', data=np.asarray([R_BH_unpin_0[i]]))
		snap_i.create_dataset('Coordinates_BH_unpin_1', data=np.asarray([R_BH_unpin_1[i]]))
		snap_i.create_dataset('Coordinates_BH_unpin_1_5', data=np.asarray([R_BH_unpin_1_5[i]]))
		snap_i.create_dataset('Velocities_BH_unpin_0', data=np.asarray([V_BH_unpin_0[i]]))
		snap_i.create_dataset('Velocities_BH_unpin_1', data=np.asarray([V_BH_unpin_1[i]]))
		snap_i.create_dataset('Velocities_BH_unpin_1_5', data=np.asarray([V_BH_unpin_1_5[i]]))
	

tf = time.time()

print(f'Values saved to file: values_BH_DM_{N_DM:.0e}.h5')
print('Runtime =', tf-ti)

