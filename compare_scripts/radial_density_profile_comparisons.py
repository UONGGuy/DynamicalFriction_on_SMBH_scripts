from load_modules import *

#### FILE PATHS #####

DM_1e4_pot_1e4 = get_snap_filename('../DM_halo_1e4/pot_1e4/isolated_example_setup_1e4/Output', 100)
DM_1e4_pot_1e5 = get_snap_filename('../DM_halo_1e4/pot_1e5/isolated_example_setup_1e4/Output', 100)
DM_1e4_pot_1e6 = get_snap_filename('../DM_halo_1e4/pot_1e6/isolated_example_setup_1e4/Output', 100)
DM_1e5_pot_1e5 = get_snap_filename('../DM_halo_1e5/pot_1e5/isolated_example_setup_1e5/Output', 100)
DM_1e5_pot_1e6 = get_snap_filename('../DM_halo_1e5/pot_1e6/isolated_example_setup_1e5_pot_1e6/Output', 100)
DM_1e6_pot_1e6 = get_snap_filename('../DM_halo_1e6/isolated_example_setup_1e6/Output', 100)
DM_1e7_pot_1e7 = get_snap_filename('../DM_halo_1e7/isolated_example_setup_1e7/Output', 100)

#### PERCENTAGE DIFFERENCE PLOTS #####

######## DM_1e4 ########

bin_1e4 = 0.5
soft_length_1e4 = 0.95
M_1e4 = np.full(int(1e4), M200/1e4)
R_1e4_pot_1e4 = np.asarray(get_snap_data(DM_1e4_pot_1e4, 1, 'Coordinates'))
R_1e4_pot_1e5 = np.asarray(get_snap_data(DM_1e4_pot_1e5, 1, 'Coordinates'))
R_1e4_pot_1e6 = np.asarray(get_snap_data(DM_1e4_pot_1e6, 1, 'Coordinates'))
D_1e4_pot_1e4 = dist_CoM(find_CoM(R_1e4_pot_1e4, M_1e4), R_1e4_pot_1e4)
D_1e4_pot_1e5 = dist_CoM(find_CoM(R_1e4_pot_1e5, M_1e4), R_1e4_pot_1e5)
D_1e4_pot_1e6 = dist_CoM(find_CoM(R_1e4_pot_1e6, M_1e4), R_1e4_pot_1e6)
rho_1e4_pot_1e4 = radial_density_halo(D_1e4_pot_1e4[:,3], M_1e4, bin_1e4, R200)
rho_1e4_pot_1e5 = radial_density_halo(D_1e4_pot_1e5[:,3], M_1e4, bin_1e4, R200)
rho_1e4_pot_1e6 = radial_density_halo(D_1e4_pot_1e6[:,3], M_1e4, bin_1e4, R200)
HQ_1e4 = hq.density(rho_1e4_pot_1e4[:,0])
diff_1e4_pot_1e4 = np.subtract(rho_1e4_pot_1e4[:,1], HQ_1e4)
diff_1e4_pot_1e5 = np.subtract(rho_1e4_pot_1e5[:,1], HQ_1e4)
diff_1e4_pot_1e6 = np.subtract(rho_1e4_pot_1e6[:,1], HQ_1e4)
cf_1e4_pot_1e4 = 100 * diff_1e4_pot_1e4 / HQ_1e4
cf_1e4_pot_1e5 = 100 * diff_1e4_pot_1e5 / HQ_1e4
cf_1e4_pot_1e6 = 100 * diff_1e4_pot_1e6 / HQ_1e4

fig_1e4 = plt.figure()
plt.plot(np.log10(rho_1e4_pot_1e4[:,0]), np.log10(HQ_1e4), 'k', label='Hernquist')
plt.plot(np.log10(rho_1e4_pot_1e4[:,0]), np.log10(rho_1e4_pot_1e4[:,1]), label='Pot. Part = 1e4')
plt.plot(np.log10(rho_1e4_pot_1e5[:,0]), np.log10(rho_1e4_pot_1e5[:,1]), label='Pot. Part = 1e5')
plt.plot(np.log10(rho_1e4_pot_1e6[:,0]), np.log10(rho_1e4_pot_1e6[:,1]), label='Pot. Part = 1e6')
plt.title('Plot of $\log_{10}$ binned halo radial density vs Hernquist profile \n(snap_100, DM Part.=1e4 [IC])')
plt.xlabel(r'$\log_{10}(r\; /[kpc])$')
plt.ylabel(r'$\log_{10}(\rho(r)\; /[10^{10}M_\odot kpc^{-3}])$')
plt.vlines(x=np.log10(2.8 * soft_length_1e4), ymin=-7.5, ymax=0, color='darkorchid', ls='dashed', label='Numerical heating limit')
plt.vlines(x=np.log10(R200), ymin=-7.5, ymax=0, color='darkblue', ls='dashed', label='R200')
plt.legend(loc='upper right')
fig_1e4.savefig('./IC_DensityProfile_DM_1e+04_pot_all_snap_100.pdf', dpi=600, bbox_inches='tight')

fig_1e4_cf = plt.figure()
plt.title('Percentage difference between simulated and Hernquist density profiles \n(snap_100 DM Part.=1e4 [IC])')
plt.xlabel(r'$\log_{10}(r\; /(kpc))$')
plt.ylabel(r'$\frac{\rho_{sim}-\rho_{H}}{\rho_{H}}\times100\%$')
plt.plot(np.log10(rho_1e4_pot_1e4[:,0]), cf_1e4_pot_1e4, label='Pot. Part = 1e4')
plt.plot(np.log10(rho_1e4_pot_1e5[:,0]), cf_1e4_pot_1e5, label='Pot. Part = 1e5')
plt.plot(np.log10(rho_1e4_pot_1e6[:,0]), cf_1e4_pot_1e6, label='Pot. Part = 1e6')
plt.vlines(x=np.log10(2.8 * soft_length_1e4), ymin=min(cf_1e4_pot_1e4), ymax=max(cf_1e4_pot_1e4), color='darkorchid', ls='dashed', label='Numerical heating limit')
plt.vlines(x=np.log10(R200), ymin=min(cf_1e4_pot_1e4), ymax=max(cf_1e4_pot_1e4), color='darkblue', ls='dashed', label='R200')
plt.legend(loc='lower left')
plt.axhline(0, 0, R200 * 1.1, ls='--', c='k')
plt.xlim([np.log10(2.8 * soft_length_1e4 * 0.9), np.log10(R200 * 1.1)])

fig_1e4_cf.savefig('./IC_DensityCompare_DM_1e+04_pot_all_snap_100.pdf', dpi=600, bbox_inches='tight')

######## DM_1e5 ########

bin_1e5 = 0.5 
soft_length_1e5 = 0.44
M_1e5 = np.full(int(1e5), M200/1e5)
R_1e5_pot_1e5 = np.asarray(get_snap_data(DM_1e5_pot_1e5, 1, 'Coordinates'))
R_1e5_pot_1e6 = np.asarray(get_snap_data(DM_1e5_pot_1e6, 1, 'Coordinates'))
D_1e5_pot_1e5 = dist_CoM(find_CoM(R_1e5_pot_1e5, M_1e5), R_1e5_pot_1e5)
D_1e5_pot_1e6 = dist_CoM(find_CoM(R_1e5_pot_1e6, M_1e5), R_1e5_pot_1e6)
rho_1e5_pot_1e5 = radial_density_halo(D_1e5_pot_1e5[:,3], M_1e5, bin_1e5, R200)
rho_1e5_pot_1e6 = radial_density_halo(D_1e5_pot_1e6[:,3], M_1e5, bin_1e5, R200)
HQ_1e5 = hq.density(rho_1e5_pot_1e5[:,0])
diff_1e5_pot_1e5 = np.subtract(rho_1e5_pot_1e5[:,1], HQ_1e5)
diff_1e5_pot_1e6 = np.subtract(rho_1e5_pot_1e6[:,1], HQ_1e5)
cf_1e5_pot_1e5 = 100 * diff_1e5_pot_1e5 / HQ_1e5
cf_1e5_pot_1e6 = 100 * diff_1e5_pot_1e6 / HQ_1e5

fig_1e5 = plt.figure()
plt.plot(np.log10(rho_1e5_pot_1e5[:,0]), np.log10(HQ_1e5), 'k', label='Hernquist')
plt.plot(np.log10(rho_1e5_pot_1e5[:,0]), np.log10(rho_1e5_pot_1e5[:,1]), label='Pot. Part = 1e5')
plt.plot(np.log10(rho_1e5_pot_1e6[:,0]), np.log10(rho_1e5_pot_1e6[:,1]), label='Pot. Part = 1e6')
plt.title('Plot of $\log_{10}$ binned halo radial density vs Hernquist profile \n(snap_100, DM Part.=1e5 [IC])')
plt.xlabel(r'$\log_{10}(r\; /(kpc))$')
plt.ylabel(r'$\log_{10}(\rho(r)\; /[10^{10}M_\odot kpc^{-3}])$')
plt.vlines(x=np.log10(2.8 * soft_length_1e5), ymin=-7.5, ymax=0, color='darkorchid', ls='dashed', label='Numerical Heating Limit')
plt.vlines(x=np.log10(R200), ymin=-7.5, ymax=0, color='darkblue', ls='dashed', label='R200')
plt.legend(loc='upper right')

fig_1e5.savefig('./IC_DensityProfile_DM_1e+05_pot_all_snap_100.pdf', dpi=600, bbox_inches='tight')

fig_1e5_cf = plt.figure()
plt.title('Percentage difference between simulated and Hernquist density profiles \n(snap_100 DM Part.=1e5 [IC])')
plt.xlabel(r'$\log_{10}(r\; /[kpc])$')
plt.ylabel(r'$\frac{\rho_{sim}-\rho_{H}}{\rho_{H}}\times100\%$')
plt.plot(np.log10(rho_1e5_pot_1e5[:,0]), cf_1e5_pot_1e5, label='Pot. Part = 1e5')
plt.plot(np.log10(rho_1e5_pot_1e6[:,0]), cf_1e5_pot_1e6, label='Pot. Part = 1e6')
plt.vlines(x=np.log10(2.8 * soft_length_1e5), ymin=min(cf_1e5_pot_1e6), ymax=max(cf_1e5_pot_1e6), color='darkorchid', ls='dashed', label='Numerical heating limit')
plt.vlines(x=np.log10(R200), ymin=min(cf_1e5_pot_1e6), ymax=max(cf_1e5_pot_1e6), color='darkblue', ls='dashed', label='R200')
plt.legend(loc='lower left')
plt.axhline(0, 0, R200*1.1, ls='--', c='k')
plt.xlim([np.log10(2.8 * soft_length_1e5 * 0.9), np.log10(R200 * 1.1)])

fig_1e5_cf.savefig('./IC_DensityCompare_DM_1e+05_pot_all_snap_100.pdf', dpi=600, bbox_inches='tight')


plt.show()
