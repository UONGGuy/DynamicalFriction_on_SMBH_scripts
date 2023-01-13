from load_modules import *

run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

a = () #time of snapshots

#### CODE: X-Y DENSITY PROJECTION EVOLUTION IN BOX ABOUT BH ####

fig1, axes1 = plt.subplots(1,5)
fig2, axes2 = plt.subplots(1,5)
box_bin = 11
box_half_width = 2 * soft_length_BH
for i in range(5):
	try:
		fname = get_snap_filename('../Output', i*20)
		R_DM = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
		R_BH = np.asarray(get_snap_data(fname, 5, 'Coordinates'))
		V_BH = np.asarray(get_snap_data(fname, 5, 'Velocities'))
		a = np.append(a, get_attribute(fname, 'Time'))
		n = get_attribute(fname, 'NumPart_ThisFile')[1]
		M = np.full(n, M200/N_DM)
		x_bin = np.linspace(R_BH[0][0] - box_half_width, R_BH[0][0] + box_half_width, box_bin)
		y_bin = np.linspace(R_BH[0][1] - box_half_width, R_BH[0][1] + box_half_width, box_bin)
		R_DM_xy_ii = np.where(np.abs(R_DM[:,2] - R_BH[0][2]) < box_half_width)
		R_DM_xy = R_DM[R_DM_xy_ii[0]]
		xy_dens, xedge, yedge, binno = stats.binned_statistic_2d(R_DM_xy[:,0], R_DM_xy[:,1], M[:len(R_DM_xy)], statistic='sum', bins=[x_bin, y_bin])
		mean_dens = np.sum(xy_dens) / (2 * box_half_width)**2
		xy_dens_rel = np.divide(xy_dens, mean_dens)
		XX, YY = np.meshgrid(x_bin, y_bin)#(xedge, yedge)
		scale_factor = 0.5 * soft_length_BH / max(np.abs(V_BH[0][0]), np.abs(V_BH[0][1]))
		dens_raw = axes1[i].pcolormesh(XX, YY, xy_dens_rel.T)
		cbar_raw = fig1.colorbar(dens_raw, ax=axes1[i], fraction=0.046, pad=0.04)
		cbar_raw.set_label(r'$\Sigma/\overline{\Sigma}$')#(10^{10}M_{\odot}kpc^{-2})$')
		arrow_raw = axes1[i].arrow(x=R_BH[0][0], y=R_BH[0][1], dx=V_BH[0][0] * scale_factor, dy=V_BH[0][1] * scale_factor, width=0.008, head_width=4 * 0.008, color='r')
		axes1[i].axes.set_aspect('equal')
		axes1[i].set_title('t=%2.2f Gyr' %(a[i]))
		axes1[i].set_xlabel('x/(kpc)')
		axes1[i].set_ylabel('y/(kpc)')
		xy_dens_filter = gaussian_filter(xy_dens, 2 * box_half_width * 1.0)
		mean_dens_filter = np.sum(xy_dens_filter) / (2 * box_half_width)**2
		xy_cf = np.divide(xy_dens_filter,mean_dens_filter)
		dens_cf = axes2[i].pcolormesh(XX, YY, xy_cf.T)
		cbar_cf = fig2.colorbar(dens_cf, ax=axes2[i], fraction=0.046, pad=0.04)
		cbar_cf.set_label(r'$\Sigma(\theta)/\overline{\Sigma}(\theta)$')
		arrow_cf = axes2[i].arrow(x=R_BH[0][0], y=R_BH[0][1], dx=V_BH[0][0] * scale_factor, dy=V_BH[0][1] * scale_factor, width=0.008, head_width=4 * 0.008, color='r')
		axes2[i].axes.set_aspect('equal')
		axes2[i].set_title('t=%2.2f Gyr' %(a[i]))
		axes2[i].set_xlabel('x/(kpc)')
		axes2[i].set_ylabel('y/(kpc)')

	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break

fig1.suptitle('SURFACE DENSITY VARIATION')
fig2.suptitle('SURFACE DENSITY VARIATION (SMOOTHED)')

fig1.set_size_inches(26, 6)
fig2.set_size_inches(26, 6)

fig1.tight_layout()
fig2.tight_layout()

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./BH_plots/BH_wake_evolution_DM_' + f'{N_DM:.0e}' + '_evolution.pdf', dpi=600, bbox_inches='tight')
	fig2.savefig('./BH_plots/BH_wake_evolution_norm_DM_' + f'{N_DM:.0e}' + '_evolution.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
	plt.show()
