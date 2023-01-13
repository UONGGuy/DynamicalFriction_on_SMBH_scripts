from load_modules import *
from config import *

#run_type = 2 #0 = save only / 1 = show only / 2 = save and show

#### HOLDERS ####

a = () # time of snapshots

#### CODE: EVOLUTION OF RADIAL DENSITY OVER TIME ####

fig = plt.figure()
radius = np.arange(0, R200 * 1.1, bin_width)
plt.plot(np.log10(radius), np.log10(hq.density(radius)), 'k', label='Hernquist')
#plt.title(r'Time evolution of $\log_{10}$ binned halo radial density vs'+'\nHernquist profile for halo (DM Part=%.0e, Pot. Part=%.0e) [IC]' %(N_DM, N_PotPart))
plt.xlabel(r'$\log_{10}(r\; [kpc])$', fontsize=18)
plt.ylabel(r'$\log_{10}(\rho(r)\; [10^{10}M_\odot /kpc^{3}])$', fontsize=18)

for i in range(11):
        fname = get_snap_filename('../Output', i*10)
        R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
        n = get_attribute(fname, 'NumPart_ThisFile')[1]
        M = np.full(n, M200/N_DM) #indiv DM particle mass assuming all identical
        R_CoM = find_CoM(R, M)
        a = np.append(a, get_attribute(fname, 'Time'))
        D = dist_CoM(R_CoM, R)
        r_density_R200 = radial_density_halo(D[:,3], M, bin_width, R200 * 1.1)
        plt.plot(np.log10(r_density_R200[:,0]), np.log10(r_density_R200[:,1]), '--', label='t=%2.2f Gyr' %(a[i]))

plt.axvline(x=np.log10(2.8 * soft_length), ls='--', color='darkorchid', label='Numerical heating limit')
plt.legend(loc='lower left', fontsize=14)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.xlim([-2.1, 1.7])
#plt.ylim([-7.25, 0])

fig.set_size_inches(14, 10)

#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
       fig.savefig('./IC_halo_plots/IC_DensityEvolution_DM_' + f'{N_DM:.0e}' + '.pdf', dpi=600, bbox_inches='tight')

if run_type == 1 or run_type == 2:
        plt.show()

