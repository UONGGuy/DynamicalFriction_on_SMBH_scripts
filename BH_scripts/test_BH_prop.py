from load_modules import *

snap_dir ="../arepo/"
snap_number = 71 #choose final snap = 100
seed_mass = 1.0e8 #choose M_BH = 1e8
snap_name = "%ssnap_%03d.hdf5" %(snap_dir,snap_number)

R_DM = np.asarray(get_snap_data(snap_name, 1, 'Coordinates'))
n = get_attribute(snap_name, 'NumPart_ThisFile')[1]
M = np.full(n, M200 /N_DM)
R_DM_CoM = find_CoM(R_DM, M)
pos5 = R_init + R_DM_CoM
print('R_DM_CoM =', R_DM_CoM)
print('pos5 =', pos5)

R_BH = np.asarray(get_snap_data(snap_name, 5, 'Coordinates'))
print('R_BH =', R_BH)

V_DM = np.asarray(get_snap_data(snap_name, 1, 'Velocities'))
D = dist_CoM(R_DM_CoM, R_DM)
V_DM_shell_ii = np.where((D[:,3] > np.linalg.norm(R_init)-0.5) & (D[:,3] < np.linalg.norm(R_init)+0.5))
V_DM_shell_mean = np.mean(np.linalg.norm(V_DM[V_DM_shell_ii], axis=1))
vel5 = np.asarray([[0, V_DM_shell_mean, 0]], dtype='float32')
print('vel5 =', vel5)

V_BH = np.asarray(get_snap_data(snap_name, 5, 'Velocities'))
print('V_BH =', V_BH)
