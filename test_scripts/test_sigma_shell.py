from load_modules import *
from config import *

snap_dir = "../arepo/"
snap_number = 100
snap_name = "%ssnap_%03d.hdf5" %(snap_dir,snap_number)
print(snap_name)

R_DM = np.asarray(get_snap_data(snap_name, 1, 'Coordinates'))
V_DM = np.asarray(get_snap_data(snap_name, 1, 'Velocities'))
n = get_attribute(snap_name, 'NumPart_ThisFile')[1]
M = np.full(n, M200 / N_DM)
R_DM_CoM = find_CoM(R_DM, M)
D = dist_CoM(R_DM_CoM, R_DM)
V_DM_mag = np.linalg.norm(V_DM, axis=1)
V_DM_shell_ii = np.where((D[:,3] > np.linalg.norm(R_init)-0.5) & (D[:,3] < np.linalg.norm(R_init)+0.5))
sigma_3 = np.std(V_DM[V_DM_shell_ii], axis=0)
sigma_r = np.average(sigma_3)
print(sigma_r)
