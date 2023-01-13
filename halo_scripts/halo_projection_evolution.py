from load_modules import *

#run_type = 2 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

a = () #time of snapshots
n_bins = 78
box_bin = np.linspace(-1.1 * R200, 1.1 * R200, n_bins)

#### CODE: X-Y DENSITY PROJECTION EVOLUTION ABOUT VIRIAL HALO ####

fig, axes = plt.subplots(1,3)#, sharex=True, sharey=True)

for i, ax in enumerate(axes.flat):
	try:
		fname = get_snap_filename('../Output', i*50)
		R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
		n = get_attribute(fname, 'NumPart_ThisFile')[1]
		a = np.append(a, get_attribute(fname, 'Time'))
		M = np.full(n, M200/N_DM)
		R_CoM = find_CoM(R, M)
		D_CoM = dist_CoM(R_CoM, R)
		M_core, bin_edge, freq = stats.binned_statistic(D_CoM[:,3], M, 'sum', bins=[0, 0.1 * R200])
		rho_core = M_core[0] / (4/3 * np.pi * (0.1 * R200) **3)
		xy_dens, xedge, yedge, binno = stats.binned_statistic_2d(D_CoM[:,0], D_CoM[:,1], M, statistic='sum', bins = [box_bin, box_bin])
		XX, YY = np.meshgrid(xedge, yedge)
		im = axes[i].pcolormesh(XX, YY, xy_dens.T)
		axes[i].axes.set_aspect('equal')
		axes[i].set_title('t=%2.2f Gyr' %(a[i]))
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break

fig.subplots_adjust(right=0.8)
cbar = fig.add_axes([0.85, 0.15, 0.015, 0.7])
fig.colorbar(im, cax=cbar, fraction=0.046, pad=0.04, label=r'$\Sigma\; [10^{10}M_{\odot}\,kpc^{-2}$]', norm=colors.LogNorm())

#fig.suptitle('Evolution of halo projected mass with time in the x-y plane (DM Part.=%.0e, Pot. Part.=%.0e) [IC]' %(N_DM, N_PotPart))
fig.supxlabel(r'$x\; [kpc]$')
fig.supylabel(r'$y\; [kpc]$')

fig.set_size_inches(14, 5)

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig.savefig('./IC_halo_plots/IC_DensityProjection_DM_' + f'{N_DM:.0e}' + '_evolution.pdf', dpi=600, bbox_inches='tight')
        
if run_type == 1 or run_type == 2:
        plt.show()
