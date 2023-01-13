from load_modules import *
from config import *

#run_type = 2 #0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

#i = 100 # snap of interest
fname = get_snap_filename('../Output', snap_i)
a = get_attribute(fname,'Time')
snap_no = f'{snap_i:03d}'
R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
n = get_attribute(fname,'NumPart_ThisFile')[1]
M = np.full(n,M200 /N_DM)
R_CoM = find_CoM(R, M)
D = dist_CoM(R_CoM, R)

#### CODE: RADIAL DENSITY PLOTS ####

######## LOG-LOG PLOT ########

radius = np.arange(0, R200 * 1.1, bin_width)
HQ_profile = hq.density(radius)

fig4 = plt.figure()
plt.plot(np.log10(radius), np.log10(HQ_profile), 'k')
r_density_R200 = radial_density_halo(D[:,3], M, bin_width, R200)
plt.plot(np.log10(r_density_R200[:,0]),np.log10(r_density_R200[:,1]),'r')

#plt.title('Plot of $\log_{10}$ binned halo radial density vs Hernquist profile \n(snap_' + snap_no + ', DM Part.= %.0e, Pot. Part=%.0e [IC])' %(N_DM, N_PotPart))
plt.xlabel(r'$\log_{10}(r\; [kpc])$', fontsize=16)
plt.ylabel(r'$\log_{10}(\rho(r)\; [10^{10}M_\odot / kpc^{3}])$', fontsize=16)
plt.vlines(x=np.log10(2.8 * soft_length), ymin=np.log10(hq.density(R200 * 1.5)), ymax=np.log10(hq.density(bin_width / 2)), color='darkorchid', ls='dotted')
plt.legend(['Hernquist', 'R200', 'Numerical heating limit'], loc='upper right')

######## PERCENTAGE DIFFERENCE PLOTS ########

r_density_R200_20 = radial_density_halo(D[:,3], M, 0.25, R200)
HQ_profile_2 = hq.density(r_density_R200_20[:,0])

density_diff = np.subtract(r_density_R200_20[:,1], HQ_profile_2)
density_cf = [a / b * 100 for a,b in zip(density_diff, HQ_profile_2)]
density_cf = density_diff / HQ_profile_2
density_cf = [x * 100 for x in density_cf]

fig5 = plt.figure()
#plt.title('Percentage difference between simulated and Hernquist density profiles \n(snap_' + snap_no + ', DM Part.=%.0e, Pot. Part=%.0e [IC])' %(N_DM, N_PotPart))
plt.xlabel(r'$\log_{10}(r\; [kpc])$', fontsize=16)
plt.ylabel(r'$(\rho_{sim}-\rho_{H})/{\rho_{H}}\times100\%$', fontsize=16)
plt.plot(np.log10(r_density_R200_20[:,0]), density_cf)
plt.vlines(x=np.log10(2.8 * soft_length), ymin=min(density_cf), ymax=max(density_cf), color='darkorchid', ls='dashed')
plt.vlines(x=np.log10(R200), ymin=min(density_cf), ymax=max(density_cf), color='darkblue', ls='dashed')
plt.xlim([np.log10(2.8 * soft_length * 0.9), np.log10(R200 * 1.1)])
plt.legend([r'$(\rho_{sim}-\rho_{H})/{\rho_{H}}\times100\%$', 'Numerical heating limit'], loc='upper right')
plt.axhline(0, 0, R200 * 1.1, ls='--', c='k')

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig4.savefig('./IC_halo_plots/IC_DensityProfile_DM_' + f'{N_DM:.0e}' + '_snap_' + snap_no + '.pdf', dpi=600)
	fig5.savefig('./IC_halo_plots/IC_DensityDiff_DM_' + f'{N_DM:.0e}' + '_snap_' + snap_no + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()

