import os
from distutils.dir_util import copy_tree
from plot_paths import *

Runs = [DM_1e4, DM_1e5, DM_1e6, DM_1e7, DM_1e8]

if not os.path.exists('../Final_Plots_Report'):
	for k in range(len(Runs)-1):
		os.makedirs(f'../Final_Plots_Report/DM_{Runs[k].N_DM:.0e}/halo_plots/')
		os.makedirs(f'../Final_Plots_Report/DM_{Runs[k].N_DM:.0e}/BH_plots/')
	os.makedirs(f'../Final_Plots_Report/DM_compare/')
	
for i in range(len(Runs)):
	copy_tree(Runs[i].halo_path, f'../Final_Plots_Report/DM_{Runs[i].N_DM:.0e}/halo_plots/')
	copy_tree(Runs[i].BH_path, f'../Final_Plots_Report/DM_{Runs[i].N_DM:.0e}/BH_plots/')

copy_tree('../../DM_halo_compare/BH_plots_compare/', '../Final_Plots_Report/DM_compare/')
