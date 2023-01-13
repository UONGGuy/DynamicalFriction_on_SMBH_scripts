from load_modules import *

#run_type = 1 # indiv script selection 0 = save only / 1 = show only / 2 = save and show

CoM_trace = np.empty([1,4])# np.empty(4) # first 3 columns = coords, 4th = time/scale factor

i = 0
while True:
	try:
		fname = get_snap_filename('../Output', i*10)
		R = np.asarray(get_snap_data(fname, 1, 'Coordinates'))
#		print(fname)
		n = get_attribute(fname,'NumPart_ThisFile')[1]
		M = np.full(n, 1)
		R_CoM = find_CoM(R, M)
		a = get_attribute(fname,'Time')
		stack = np.hstack((R_CoM,a))
		CoM_trace = np.append(CoM_trace, [stack], axis=0)
#		print(CoM_trace[i])
	except(KeyError, OSError, NameError, UnboundLocalError, IOError):
		break
	else:
                i += 1
#print(CoM_trace)

#projection in xy, xz, yz planes

fig1 = plt.figure()
ax = fig1.add_subplot(projection='3d')
ax.scatter3D(CoM_trace[:,0], CoM_trace[:,1], CoM_trace[:,2])
plt.title('Plot tracing halo CoM location through time [IC]')
ax.set_xlabel(r'$x\; [kpc]$')
ax.set_ylabel(r'$y\; [kpc]$')
ax.set_zlabel(r'$z\; [kpc]$')


#### BULK RUN OPTIONS ####

if run_type == 0 or run_type == 2:
	fig1.savefig('./IC_halo_plots/IC_halo_CoM_trace.pdf')        

if run_type == 1 or run_type == 2:
        plt.show()


