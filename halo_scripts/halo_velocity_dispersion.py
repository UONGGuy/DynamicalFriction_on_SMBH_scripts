from load_modules import *
from config import *

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

#snap_i = 100 #snap of interest
#np.linalg.norm(R_init) = # 0.2*R200 = 7 kpc # initial BH orbital radius
a = () # time of snapshots
fname = get_snap_filename('../Output/', snap_i)
a_i = get_attribute(fname,'Time')
R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
n = get_attribute(fname,'NumPart_ThisFile')[1]
M = np.full(n, M200 / N_DM)
R_CoM = find_CoM(R, M)
D = dist_CoM(R_CoM, R)
V_3D = np.asarray(get_snap_data(fname, 1, 'Velocities'))
radial_bin = np.arange(0, R200 * 1.1, 1)

#### CODE: DM HALO VELOCITY DISPERSION ####

sigma_x, x_edge, x_binno = stats.binned_statistic(D[:,3], V_3D[:,0], statistic='std', bins=radial_bin)
sigma_y, y_edge, y_binno = stats.binned_statistic(D[:,3], V_3D[:,1], statistic='std', bins=radial_bin)
sigma_z, z_edge, z_binno = stats.binned_statistic(D[:,3], V_3D[:,2], statistic='std', bins=radial_bin)
sigma_3d = np.c_[sigma_x, sigma_y, sigma_z]
sigma_mean = np.mean(sigma_3d, axis=1)

fig1 = plt.figure()

plt.bar(radial_bin[:-1], sigma_mean, width=np.diff(x_edge) * 0.86, align='edge', alpha=0.7, label=r'$\overline{\sigma}_{DM}$')
plt.plot(radial_bin, hq.v_disp_radial(radial_bin), label=r'$\sigma_r$', c='r')

#plt.title(f'Mean velocity dispersion of bins in halo \n(snap_{snap_i:.03d}, DM Part.={N_DM:.0e}, Pot. Part={N_PotPart:.0e} [IC])')
plt.xlabel(r'$r\; [kpc]$', fontsize=16)
plt.ylabel(r'$\sigma\; [km/s]$', fontsize=16)

plt.legend(loc='upper right', fontsize=16)


#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
        fig1.savefig(f'./IC_halo_plots/IC_halo_velocity_disp_DM_{N_DM:.0e}_snap_{snap_i:03d}.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
        plt.show()

