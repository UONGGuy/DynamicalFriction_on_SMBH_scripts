from load_modules import *

#### HOLDERS ####

i = 00 #snap of interest
a = () # time of snapshots
fname = get_snap_filename('../Output', i)
a_i = get_attribute(fname,'Time')
snap_no = f'{i:03d}'
R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
n = get_attribute(fname,'NumPart_ThisFile')[1]
M = np.full(n, M200 / N_DM)
R_CoM = find_CoM(R, M)
D = dist_CoM(R_CoM, R)
V_3D = np.asarray(get_snap_data(fname, 1, 'Velocities'))
V_BH = np.asarray(get_snap_data(fname, 5, 'Velocities'))

#### CODE: SORTING METHOD ####

#R_i = [7 / np.sqrt(3), 7 / np.sqrt(3), 7 / np.sqrt(3)]
D_i = D[:,:-1] - R_init
D_i_Mag = np.linalg.norm(D_i, axis=1)
R_V_i = np.c_[D_i_Mag, V_3D]
R_V_i = R_V_i.tolist()
R_V_i = sorted(R_V_i, key=lambda x: x[0])
R_V_i = np.asarray(R_V_i)
R_V_i_64 = R_V_i[0:64,:]

M_64 = M[0:64]
V_CoM = find_CoM(R_V_i_64[:,1:4], M_64)
V_rel_i = np.subtract(R_V_i_64[:,1:4], V_CoM)
R_V_rel_i_64 = np.c_[R_V_i_64[:,0], np.linalg.norm(V_rel_i, axis=1)]

#print(V_rel_i)

fig1 = plt.figure()
#plt.plot(R_V_rel_i_64[:,0], R_V_rel_i_64[:,1], 'x')
Prob, bins, patches = plt.hist(R_V_rel_i_64[:,1], 20, density=True, alpha=0.6)

def Maxwell_distrib(v, N, s):
        return N * np.sqrt(2 / np.pi) * v**2 * np.exp(-v**2 / (2 * s)) / np.sqrt(s**3)

p1, p2 = optimize.curve_fit(Maxwell_distrib, bins[:-1], Prob)

b = np.linspace(0, 80, 801)
fit1 = (p1[0] * np.sqrt(2 / np.pi) * b[:-1]**2 * np.exp(-b[:-1]**2 / (2 * p1[1])) / np.sqrt(p1[1]**3))
plt.plot(b[:-1], fit1, label='curve fit')
fit2 = (p1[0] * np.sqrt(2 / np.pi) * bins[:-1]**2 * np.exp(-bins[:-1]**2 / (2 * p1[1])) / np.sqrt(p1[1]**3))
plt.plot(bins[:-1], fit2, label='curve fit 2', c='r')

#plt.show()

#### CODE: KD-TREE METHOD ####

tree = cKDTree(D_i)
dd, ii = tree.query([0, 0, 0], k=64)
#print(dd, ii, sep='\n')

tree2 = cKDTree(D[:,:-1])
dd2, ii2 = tree2.query(R_init, k=64)
#print(dd2, ii2, sep='\n')
V_tree2 = [V_3D[index] - V_CoM for index in ii2][0]
#V_tree2 = V_tree2 - V_CoM
#V_tree2 = V_tree2[0,:]
#print('V_tree2', V_tree2)
#print(dd2[0][-1])

rho_local = (64 * M200 / N_DM) / (4 / 3 * np.pi * dd2[0][-1]**3) #in units of 10^{10}M_sun kpc^{-3}
v_local_mag = np.linalg.norm(V_tree2, axis=1)
v_mean_local = np.mean(v_local_mag) #in units of km s^{-1}
#print(rho_local)
#print(v_mean_local)
fig2 = plt.figure()
Prob2, bins2, patches2 = plt.hist(np.linalg.norm(V_tree2, axis=1), density=True)
ProbDensity = np.multiply(Prob2, np.diff(bins2))
a1, a2 = optimize.curve_fit(Maxwell_distrib, bins2[:-1], Prob2)
fit_1 = (a1[0] * np.sqrt(2 / np.pi) * b[:-1]**2 * np.exp(-b[:-1]**2 / (2 * p1[1])) / np.sqrt(p1[1]**3))
plt.plot(b[:-1], fit_1)


#plt.show()

#### CODE: LOCAL DENSITY ####

tree3 = cKDTree(D[:,:-1])
dd3, ii3 = tree3.query(R_init, k=64)
r_tree3 = D_i[ii3]
v_tree3 = V_3D[ii3] - V_BH
v_local_mag3 = np.linalg.norm(V_tree3, axis=1)
v_thresh_ii = np.where(v_local_mag3 < np.linalg.norm(V_BH))
v_thresh = v_local_mag[v_thresh_ii]
rho_local = len(ii3) * M200 / N_DM / (4 / 3 * np.pi * max(dd3[0]))
v_local_mean = np.mean(np.linalg.norm(v_tree3, axis=1))
v_thresh_mean = np.mean(v_thresh)

print('dd3\n', max(dd3[0]))
print('v_tree3\n', v_tree3)
print('v_BH\n', np.linalg.norm(V_BH))
print('v_thresh_ii\n', v_thresh_ii)
print('v_thresh\n', v_thresh)
print('rho_locali\n', rho_local)
print('v_local_mean\n', v_local_mean)
print('v_thresh_mean\n', v_thresh_mean)


