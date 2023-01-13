from load_modules import *
import h5py 

f = h5py.File(f'values_BH_DM_{N_DM:.0e}.h5', 'r')

ls1 = list(f.keys()) # List of keys in this directory
#print(ls1)

i = 0

#print(f['/snap_%d' %(i)].attrs['Time'])
#print(get_snap_attr_2('values_BH.h5', i, 'Time'))
#print(f['/snap_100'].attrs['NumPart_DM'])
R_BH = f['/snap_000/Coordinates_BH']
print(np.asarray(R_BH)[0,3])
print(get_snap_data_2(f'values_BH_DM_{N_DM:.0e}.h5', i, 'Coordinates_BH')[0,3])
