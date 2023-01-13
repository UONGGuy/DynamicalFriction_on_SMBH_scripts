from load_modules import *
from config import *

#run_type = 2 #0 = save only / 1 = show only / 2 = save and show
#snap_i = 100
a = ()
n_bin = 71

#### CODE: EXTRACT PARTICLE LOCATIONS ####

fname = get_snap_filename('../Output', snap_i)
R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
n = get_attribute(fname, 'NumPart_ThisFile')[1]
a = np.append(a, get_attribute(fname, 'Time'))
M = np.full(n, M200 / N_DM)
R_CoM = find_CoM(R, M)
D = dist_CoM(R_CoM, R)
box_bin = np.linspace(-0.1 * R200, 0.1 * R200, n_bin)
XX, YY = np.meshgrid(box_bin, box_bin)

#### CODE: DENSITY PROJECTION PLOT ####

fig, axes = plt.subplots(1,1)

M_core, edges,  bn = stats.binned_statistic(D[:,3], M, statistic='sum', bins = [0, 0.1 * R200])
xy_dens, xedge, yedge, binno = stats.binned_statistic_2d(D[:,0], D[:,1], M, statistic='sum', bins = [box_bin, box_bin])
rho_core = M_core[0] / (4/3 * np.pi * (0.1 * R200) **3)

im = axes.pcolormesh(XX, YY, xy_dens.T, norm=colors.LogNorm())#vmin=0.01 * rho_core, vmax=1.2 * rho_core))
axes.title.set_text('Time = %2.2f GYr' %(a[0]))
axes.set_aspect('equal', adjustable='box')
Sigma = (M200 * 1e10) / ((box_bin[-1] - box_bin[0]) / (n_bin - 1))**2
cbar = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.04)
cbar.set_label(r'$\log_{10}(\Sigma/[%2.0e\; M_{\odot}\,/kpc^{2}])$' %(Sigma), fontsize=16)
axes.axes.set_xlabel(r'$x\; [kpc]$', fontsize=16)
axes.axes.set_ylabel(r'$y\; [kpc]$', fontsize=16)

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig.savefig('./IC_halo_plots/IC_DensityProjection_DM_' + f'{N_DM:.0e}' + '_snap_' + f'{snap_i:03d}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type ==1 or run_type == 2:
	plt.show()
