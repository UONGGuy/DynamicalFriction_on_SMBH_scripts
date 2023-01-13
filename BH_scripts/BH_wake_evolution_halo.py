from load_modules import *

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

fname = f'values_BH_DM_{N_DM:.0e}.h5'
fig1, ax = plt.subplots(1,3)

######## CODE: WAKE PROJECTION EVOLUTION ########

for k in range(3):
	snap_no = (k + 1) * 26
	if N_DM == 1e4:
		box_bin = 15
		box_half_width = 7.5
	elif N_DM == 1e5:
		box_bin = 41
		box_half_width = 7
	elif N_DM == 1e6:
		box_bin = 101
		box_half_width = 5.5
	elif N_DM >= 1e7:
		if k == 0:
			box_bin = 121
			box_half_width = 6
		elif k == 1:
			box_bin = 81
			box_half_width = 4
		elif k == 2:
			box_bin = 61
			box_half_width = 3
#EXTRACT DATA
	R_DM_bg = np.asarray(get_snap_data_2(fname, snap_no-1, 'Coordinates_DM'))
	R_BH_prev = np.asarray(get_snap_data_2(fname, snap_no-1, 'Coordinates_BH'))
	R_BH_prev2 = np.asarray(get_snap_data_2(fname, snap_no-2, 'Coordinates_BH'))
	R_DM = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_DM'))
	R_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Coordinates_BH'))
	V_BH = np.asarray(get_snap_data_2(fname, snap_no, 'Velocities_BH'))
	a = get_snap_attr_2(fname, snap_no, 'Time')
	n = get_snap_attr_2(fname, snap_no, 'NumPart_DM')
	M = np.full(n, M200 / N_DM)
#SET-UP CANVAS
	x_bin = np.linspace(-box_half_width, box_half_width, box_bin)
	y_bin = x_bin
	xbox, ybox = x_bin, y_bin
	XX, YY = np.meshgrid(xbox, ybox)
#SURFACE DENSITY RAW
	R_DM_xy_ii = np.where(np.abs(R_DM[:,2]) < box_half_width)
	R_DM_xy = R_DM[R_DM_xy_ii[0]]
	xy_dens, xedge, yedge, binno = stats.binned_statistic_2d(R_DM_xy[:,0], R_DM_xy[:,1], M[:len(R_DM_xy)], statistic='sum', bins=[x_bin, y_bin])
#SURFACE DENSITY COMPARISON
	R_DM_xy_ii_bg = np.where(np.abs(R_DM_bg[:,2]) < box_half_width)
	R_DM_xy_bg = R_DM_bg[R_DM_xy_ii_bg[0]]
	xy_dens_bg, xedge_bg, yedge_bg, binno_bg = stats.binned_statistic_2d(R_DM_xy_bg[:,0], R_DM_xy_bg[:,1], M[:len(R_DM_xy_bg)], statistic='sum', bins=[x_bin, y_bin])
	xy_diff = xy_dens - xy_dens_bg
	xy_cf = np.divide(xy_diff, xy_dens_bg)
	xy_cf = gaussian_filter(xy_cf, sigma = box_half_width * 0.90)
	dens_cf = ax[k].pcolormesh(XX, YY, xy_cf.T)
	cbar_cf = fig1.colorbar(dens_cf, ax=ax[k], fraction=0.046, pad=0.04)
	cbar_cf.set_label(r'$(\Sigma-\Sigma_0)/\Sigma_0$', size=7.5)
	cbar_cf.ax.tick_params(labelsize=6)
	ax[k].axes.set_aspect('equal')
	ax[k].set_title(f'Time={a:2.2f} Gyr', size=8)
	ax[k].axes.set_xlabel('x [kpc]', size=7.5)
	ax[k].axes.set_ylabel('y [kpc]', size=7.5)
	ax[k].tick_params(axis='both', labelsize=6)
	BH_trajectory = ax[k].plot([R_BH_prev2[0][0], R_BH_prev[0][0], R_BH[0][0]], [R_BH_prev2[0][1], R_BH_prev[0][1], R_BH[0][1]], 'r')
	BH_loc_curr = ax[k].plot(R_BH[0][0], R_BH[0][1], '.k')

fig1.tight_layout()
#fig1.subplots_adjust(right=0.8)
#cbar_ax = fig1.add_axes([0.85, 0.35, 0.02, 0.3])
#fig1.colorbar(dens_cf, cax=cbar_ax, label=r'$(\Sigma-\Sigma_0)/\Sigma_0$')

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./BH_plots/BH_wake_DM_' + f'{N_DM:.0e}' + '_evolution.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()
