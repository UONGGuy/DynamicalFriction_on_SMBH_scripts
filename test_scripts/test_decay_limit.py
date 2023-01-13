from load_modules import *
from config import *

num_lim = 2.8 * soft_length
R_BH = np.linalg.norm(R_init)
snap_no = 0

while R_BH > num_lim:
	snap_no += 1
	fname = f'values_BH_DM_{N_DM:.0e}.h5'
	R_BH = get_snap_data_2(fname, snap_no, 'Coordinates_BH')[0,3]
	
print(snap_no)
t_decay = get_snap_attr_2(fname, snap_no, 'Time')
print(t_decay)

