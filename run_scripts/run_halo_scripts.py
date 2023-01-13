import os
import time

ti = time.time()

if not os.path.exists('IC_halo_plots'):
    os.makedirs('IC_halo_plots')

import halo_CoM_trace
import halo_density_single
import halo_density_evolution
import halo_part_distrib
import halo_projection_single
import halo_projection_evolution
import halo_velocity_dispersion 

tf = time.time()

print('\n')
print(f'Time for run = {tf-ti:.2f}s')
print('See config.py for settings and values.')
print('\n')

#with open("./config.py", "r") as f:
#	print(f.read())


