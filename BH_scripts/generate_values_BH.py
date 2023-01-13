from load_modules import *

#runtime ~9mins

ti = time.time()

hf = h5py.File('values_BH_DM_' + f'{N_DM:.0e}' + '.h5', 'w')
i = 0

while True:
	try:
		fname = get_snap_filename('../Output', i)
		snap_i = hf.create_group(f'snap_{i:03d}')
		a = get_attribute(fname, 'Time')
		snap_i.attrs['Time'] = a
		n = get_attribute(fname,'NumPart_ThisFile')[1]
		snap_i.attrs['NumPart_DM'] = n
		M = np.full(n, M200 /N_DM)
		R_DM = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
		R_CoM = find_CoM(R_DM, M)
		D = dist_CoM(R_CoM, R_DM)
		snap_i.create_dataset('Coordinates_DM', data=D)
		V_DM = np.asarray(get_snap_data(fname, 1, 'Velocities'))
		snap_i.create_dataset('Velocities_DM', data=V_DM)
		R_BH = np.asarray(get_snap_data(fname, 5, 'Coordinates'))
		R_BH_CoM = np.asarray(R_BH - R_CoM)
		R_BH_CoM_i = np.c_[R_BH_CoM, np.linalg.norm(R_BH_CoM)]
		snap_i.create_dataset('Coordinates_BH', data=R_BH_CoM_i)
		V_BH = np.asarray(get_snap_data(fname, 5, 'Velocities'))
		snap_i.create_dataset('Velocities_BH', data=V_BH)
		print(i)
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		tf = time.time()
		print(f'Error occured during snap_{i:3d}')
		print('Runtime = ', tf-ti)
		break
	else:
		i += 1

tf = time.time()

print('Runtime = ', tf-ti)	
