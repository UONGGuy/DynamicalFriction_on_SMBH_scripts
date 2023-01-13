from load_modules import * 
from config import *

#run_type = 2 #0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

i = 00 # snap of interest
fname = get_snap_filename('../Output', i)
a = get_attribute(fname,'Time')
snap_no = f'{i:03d}'
R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
n = get_attribute(fname,'NumPart_ThisFile')[1]
M = np.full(n, M200 / N_DM)
R_CoM = find_CoM(R, M)
D = dist_CoM(R_CoM, R)

#### CODE: PARTICLE DISTRIBUTION HISTOGRAMS ####

######## FULL DISTRIBUTION ########

fig1 = plt.figure()
part_bin_1, edges_1, bin_no_1 = stats.binned_statistic(D[:,3], M * N_DM, 'count', bins=20)
plt.bar(edges_1[:-1], part_bin_1, width=np.diff(edges_1)*0.85, align='edge', ec='k')
plt.title('Histogram of all particles in halo \n(snap_' + snap_no + ', DM Part.=%.0e, Pot. Part=%.0e [IC])' %(N_DM, N_PotPart))
plt.xlabel(r'$r\; [kpc]$')
plt.ylabel(r'Particle frequency')

######## INNER 10% OF R200 ########

fig2 = plt.figure()
part_bin_2, edges_2, bin_no_2 = stats.binned_statistic(D[:,3], M * N_DM, 'count', bins=np.linspace(0, 0.1*R200, 125))
plt.bar(edges_2[:-1], part_bin_2, width=np.diff(edges_2)*0.85, align='edge', ec='k')
plt.title('Histogram of particles within 10% of R200=35kpc \n(snap_' + snap_no + ', DM Part.=%.0e, Pot. Part=%.0e [IC])' %(N_DM, N_PotPart))
plt.xlabel(r'$r\; [kpc]$')
plt.ylabel(r'Particle frequency')

######## BIN WIDTH ESTIMATION (~20 particles enclosed) ########

fig3 = plt.figure()
part_bin_3, edges_3, bin_no_3 = stats.binned_statistic(D[:,3], M * N_DM, 'count', bins=np.linspace(0, bin_width, 125))
plt.bar(edges_3[:-1], part_bin_3, width=np.diff(edges_3)*0.85, align='edge', ec='k')
plt.title('Histogram of particles in bin width=%4.4f kpc \n(snap_'%(bin_width) + snap_no + ', DM Part.=%.0e, Pot. Part=%.0e [IC])' %(N_DM, N_PotPart))
plt.xlabel(r'$r\; [kpc]$')
plt.ylabel(r'Particle frequency')

######## SHELL THICKNESS (~500 particles enclosed) AT R_init ########

fig4 = plt.figure()
part_bin_4, edges_4, bin_no_4 = stats.binned_statistic(D[:,3], M * N_DM, 'count', bins=np.linspace(5.5, 10.5, 6))
plt.bar(edges_4[:-1], part_bin_4, width=np.diff(edges_4)*0.85, align='edge', ec='k')
R_init_mag = np.linalg.norm(R_init)
plt.title(f'Histogram of particles about {R_init_mag:.2f} kpc of halo \n(snap_{snap_no}, DM Part={N_DM:.0e}, Pot. Part={N_PotPart:.0e} ')#'Histogram of particles about %.2f kpc of halo \n(snap_' + snap_no + ', DM Part.=%.0e, Pot. Part=%.0e [IC])' %(R_init_mag, N_DM, N_PotPart))
plt.xlabel(r'$r\; [kpc]$')
plt.ylabel(r'Particle frequency')

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./IC_halo_plots/IC_PartHist_full_DM_' + f'{N_DM:.0e}' + '_snap_' + snap_no + '.pdf', dpi=600, bbox_inches='tight')
	fig2.savefig('./IC_halo_plots/IC_PartHist_0.1R200_DM_' + f'{N_DM:.0e}' + '_snap_' + snap_no + '.pdf', dpi=600)
	fig3.savefig('./IC_halo_plots/IC_PartHist_bin_width_DM_' + f'{N_DM:.0e}' + '_snap_' + snap_no + '.pdf', dpi=600)
	fig4.savefig('./IC_halo_plots/IC_PartHist_7kpc_shell_DM_' + f'{N_DM:.0e}' + '_snap_' + snap_no + '.pdf', dpi=600)

if run_type == 1 or run_type == 2:
	plt.show()
