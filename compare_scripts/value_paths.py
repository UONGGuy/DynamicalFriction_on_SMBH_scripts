class Values:
	def __init__(self, N_DM, halo_path, BH_path, soft_length):
		self.N_DM = N_DM
		self.halo_path = halo_path
		self.BH_path = BH_path
		self.soft_length = soft_length

DM_1e4 = Values(1e4, '../DM_halo_1e4/IC_halo/isolated_example_setup_1e4/Output/', '../DM_halo_1e4/IC_BH/isolated_example_setup_1e4/python_scripts/values_BH_DM_1e+04.h5', 0.95)
DM_1e5 = Values(1e5, '../DM_halo_1e5/IC_halo/isolated_example_setup_1e5/Output/', '../DM_halo_1e5/IC_BH/isolated_example_setup_1e5/python_scripts/values_BH_DM_1e+05.h5', 0.44)
DM_1e6 = Values(1e6, '../DM_halo_1e6/IC_halo/isolated_example_setup_1e6/Output/', '../DM_halo_1e6/IC_BH/isolated_example_setup_1e6/python_scripts/values_BH_DM_1e+06.h5', 0.20)
DM_1e7 = Values(1e7, '../DM_halo_1e7/IC_halo/isolated_example_setup_1e7/Output/', '../DM_halo_1e7/IC_BH/isolated_example_setup_1e7/python_scripts/values_BH_DM_1e+07.h5', 0.09)
DM_1e8 = Values(1e8, '../DM_halo_1e8/IC_halo/isolated_example_setup_1e8/Output/', '../DM_halo_1e8/IC_BH/isolated_example_setup_1e8/python_scripts/values_BH_DM_1e+08.h5', 0.042)
