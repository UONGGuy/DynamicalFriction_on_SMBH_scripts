from load_modules import *
from config import *

with open('BH_analytic_decay_values.csv') as f:
	csvreader = csv.reader(f)
	header = next(csvreader)
	a = [np.float32(j) for j in next(csvreader)]
	R_BH_list = next(csvreader)
	R_BH_CoM = np.array([np.array(ast.literal_eval(R_BH_list[j])) for j in range(len(R_BH_list))])
	V_BH_list = next(csvreader)
	rho_local = np.asarray([np.float64(k) for k in next(csvreader)])

fname = f'values_BH_DM_{N_DM:.0e}.h5'
snap_no = 0

rho_local_h5, sigma_h5 = [], []
while True:
	try: 
		rho_local_h5 = np.append(rho_local_h5, get_snap_attr_2(fname, snap_no, 'rho_local_BH'))
		sigma_h5 = np.append(sigma_h5, get_snap_attr_2(fname, snap_no, 'sigma_local_BH'))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else: snap_no += 1

print('sigma', sigma_h5, sep='\n')
print(rho_local, rho_local_h5, sep='\n')
print(type(rho_local), type(rho_local_h5))

cf = np.subtract(rho_local, rho_local_h5)
print(cf)
