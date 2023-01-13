import os
import time
from config import N_DM

ti = time.time()

if not os.path.exists('BH_plots'):
    os.makedirs('BH_plots')

try:
	import BH_orbit_projection_plot
	import BH_velocity_plots
	import BH_AngMom_plots
	import BH_radial_decay_plot
	import BH_wake_evolution_halo
	if N_DM >= 1e7:
		import BH_analytic_decay_pinned_plot
		import BH_analytic_decay_unpinned_plot
except:
	print('\nPlease check if values have been generated for BH and analytic decay.\n')

tf = time.time()

print('\n')
print(f'Time for run = {tf-ti:.2f}s')
print('See config.py for settings and values.')
print('\n')

#with open("./config.py", "r") as f:
#	print(f.read())


