from load_modules import *

run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

#test snap = 76, box = 4kpc x 4kpc
fname = f'values_BH_DM_{N_DM:.0e}.h5'
box_bin = 61
snap_no = 52
box_half_width = 3
R_DM_bg = np.asarray(get_snap_data_2(fname, snap_no-1, 'Coordinates_DM'))
R_BH_prev = np.asarray(get_snap_data_2(fname, snap_no-1, 'Coordinates_BH'))
R_BH_prev2 = np.asarray(get_snap_data_2(fname, snap_no-2, 'Coordinates_BH'))
R_DM = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_DM'))
R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
a = get_snap_attr_2(fname, snap_no, 'Time')
n = get_snap_attr_2(fname, snap_no, 'NumPart_DM')
M = np.full(n, M200 / N_DM)
x_bin = np.linspace(-box_half_width, box_half_width, box_bin)
y_bin = x_bin
xbox, ybox = x_bin, y_bin
XX, YY = np.meshgrid(xbox, ybox)

#### CODE: X-Y DENSITY PROJECTIONS ####

fig1, ax = plt.subplots(1,2)

######## CODE: X-Y SURFACE DENSITY RAW ########

R_DM_xy_ii = np.where(np.abs(R_DM[:,2]) < box_half_width)
R_DM_xy = R_DM[R_DM_xy_ii[0]]
xy_dens, xedge, yedge, binno = stats.binned_statistic_2d(R_DM_xy[:,0], R_DM_xy[:,1], M[:len(R_DM_xy)], statistic='sum', bins=[x_bin, y_bin])
dens_raw = ax[0].pcolormesh(XX, YY, np.log(xy_dens).T)
cbar_raw = fig1.colorbar(dens_raw, ax=ax[0], fraction=0.046, pad=0.04)
cbar_raw.set_label(r'$\Sigma/[10^{10}M_{\odot}/kpc^{2}]$')
ax[0].axes.set_aspect('equal')
ax[0].title.set_text('SURFACE DENSITY')
ax[0].axes.set_xlabel('x/[kpc]')
ax[0].axes.set_ylabel('y/[kpc]')
BH_loc_curr = ax[0].plot(R_BH[0][0], R_BH[0][1], '.k')

######## CODE: X-Y SURFACE DENSITY COMPARISON ########

R_DM_xy_ii_bg = np.where(np.abs(R_DM_bg[:,2]) < box_half_width)
R_DM_xy_bg = R_DM_bg[R_DM_xy_ii_bg[0]]
xy_dens_bg, xedge_bg, yedge_bg, binno_bg = stats.binned_statistic_2d(R_DM_xy_bg[:,0], R_DM_xy_bg[:,1], M[:len(R_DM_xy_bg)], statistic='sum', bins=[x_bin, y_bin])
xy_diff = xy_dens - xy_dens_bg
xy_cf = np.divide(xy_diff, xy_dens_bg)
xy_cf = gaussian_filter(xy_cf, sigma = box_half_width * 0.90)
dens_cf = ax[1].pcolormesh(XX, YY, xy_cf.T)
cbar_cf = fig1.colorbar(dens_cf, ax=ax[1], fraction=0.046, pad=0.04)
cbar_cf.set_label(r'$(\Sigma-\Sigma_0)/\Sigma_0$')
ax[1].axes.set_aspect('equal')
ax[1].title.set_text('SURFACE DENSITY VARIATION')
ax[1].axes.set_xlabel('x/[kpc]')
ax[1].axes.set_ylabel('y/[kpc]')
BH_trajectory = ax[1].plot([R_BH_prev2[0][0], R_BH_prev[0][0], R_BH[0][0]], [R_BH_prev2[0][1], R_BH_prev[0][1], R_BH[0][1]], 'r')
BH_loc_curr = ax[1].plot(R_BH[0][0], R_BH[0][1], '.k')

fig1.suptitle(f'Time = {a:2.2f} Gyr (snap_{snap_no:03d})')
fig1.tight_layout()

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
        fig1.savefig('./BH_plots/BH_wake_DM_' + f'{N_DM:.0e}' + '_snap_' + f'{snap_no:03d}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
        plt.show()

