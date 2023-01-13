from load_modules import *

#### HOLDERS ####

i = 100 #snap of interest
a = () # time of snapshots
fname = get_snap_filename('../Output_no_BH', i)
a_i = get_attribute(fname,'Time')
snap_no = f'{i:03d}'
R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
n = get_attribute(fname,'NumPart_ThisFile')[1]
M = np.full(n, M200 /N_DM)
R_CoM = find_CoM(R, M)
D = dist_CoM(R_CoM, R)
V_3D = np.asarray(get_snap_data(fname, 1, 'Velocities'))
V = np.linalg.norm(V_3D, axis=1)


fig1 = plt.figure()
part_bin_1, edges_1, bin_no_1 = stats.binned_statistic(D[:,3], M * N_DM, 'count', bins=np.linspace(5.5, 10.5, 6))
plt.bar(edges_1[:-1], part_bin_1, width=np.diff(edges_1)*0.85, align='edge', ec='k')
plt.title('Histogram of all particles in halo \n(snap_' + snap_no + ', DM Part.=%.0e, Pot. Part=%.0e [IC])' %(N_DM, N_PotPart))
plt.xlabel('Distance from CoM [kpc]')
plt.ylabel('Particle frequency')

plt.show()

v_mean, edges_2, bin_no_2 = stats.binned_statistic(D[:,3], V, 'mean', bins=[6.5,7.5])

print(v_mean)
print(np.sqrt(3)*hq.v_disp_radial(7))

M_enc, edges_3, bin_no_3 = stats.binned_statistic(D[:,3], M, 'sum', bins=[0, 7])

print(M_enc)
v_calc = np.sqrt(hq.G * M_enc/7)
print(v_calc)
print(np.sqrt(hq.G * M200 * 7 / (7 + scale_length)**2))
print(hq.v_circ(7))
