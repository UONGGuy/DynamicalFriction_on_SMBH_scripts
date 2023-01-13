import numpy as np
import h5py
#import read_fof_files as rff #you won't need to load this, this is for halo catalogues from cosmological simulations
import read_snap_files as rsf
from halo_calibrate import *
#import csv
#from collate_snap_files import *
from config import *

snap_dir ="../arepo/"
snap_number = 100 #choose final snap = 100
seed_mass = 1.0e8 #choose M_BH = 1e8

#fof_name = "%sfof_subhalo_tab_%03d.hdf5" %(snap_dir,snap_number) #not needed in your case
snap_name = "%ssnap_%03d.hdf5" %(snap_dir,snap_number) 

seed_mass /= 1.0e10 #convert to 1e10 M_sun

"""Set up BH data:"""
#Read in subhalo data from cosmological simulation --> these two lines may be deleted for your example as we will set position and velocity manually (see below).
#pos_most_bound = rff.get_subhalo_data(fof_name,"SubhaloPos")
#pec_vel = rff.get_subhalo_data(fof_name,"SubhaloVel")

#Set position of BH, here I set this to the position of the most bound particle from the subhalo --> you will want to set this position manually for a suitable radius
pos5 = np.zeros((1,3),dtype="float32")
#pos5[0:,] = pos_most_bound[0]
#read in DM CoM and +7
R_DM = np.asarray(rsf.get_snap_data(snap_name, 1, 'Coordinates'))
n = rsf.get_attribute(snap_name, 'NumPart_ThisFile')[1]
M = np.full(n, M200 /N_DM)
R_DM_CoM = find_CoM(R_DM, M)
pos5 = R_init + R_CoM # 0.2 * R200 = 7 kpc

#Set velocity of BH, similarly I have set this to peculiar velocity of the halo here, while for our use case this value will be set manually.
vel5 = np.zeros(np.shape(pos5),dtype="float32")
#vel5[0:,] = pec_vel[0] #halo velocity
V_DM = np.asarray(get_snap_data(snap_name, 1, 'Velocities'))
D = dist_CoM(R_DM_CoM, R_DM)
V_DM_shell_ii = np.where((D[:,3] > np.linalg.norm(R_init)-0.5) & (D[:,3] < np.linalg.norm(R_init)+0.5))
V_DM_shell_mean = np.mean(np.linalg.norm(V_DM[V_DM_shell_ii], axis=1)) #mean magnitude of velocities
#sigma_3 = np.std(V_DM[V_DM_shell_ii], axis=0)
#sigma_r = np.average(sigma_3, weights=M[V_DM_shell_ii])
vel5 = np.asarray([[0, V_DM_shell_mean, 0]], dtype='float32')

present_part_types = [1] #we only have DM particles, in the future we will want to add to this, e.g. 0 for gas

#Loop over particle types to find maximum particle ID present
MaxIDs = []
for i in present_part_types:
    IDs= rsf.get_snap_data(snap_name,i,"ParticleIDs")
    MaxIDs.append(np.max(IDs))

#Set ID of new BH particle one higher than Max ID (ensures that this ID is not taken)
IDs5 = np.max(MaxIDs) + 1
ID5 = np.ones(len(pos5),dtype="uint32")*IDs5

partmass5 = np.ones(len(pos5),dtype="float32")*seed_mass

"""Open snap file to be modified."""
f = h5py.File(snap_name,'a')

"""Update header."""

head = f['/Header']
num_part_old = head.attrs['NumPart_ThisFile']
num_part_total_old = head.attrs['NumPart_Total']

del head.attrs['NumPart_ThisFile']
del head.attrs['NumPart_Total']

num_part_new = num_part_old
num_part_new[5] = 1

num_part_total_new = num_part_total_old
num_part_total_new[5] = 1

head.attrs.create('NumPart_ThisFile', data = num_part_new)
head.attrs.create('NumPart_Total', data = num_part_total_new)

"""Add particle group 5"""
try:
    part5 = f.create_group('PartType5')
except:
    del f['/PartType5']
    part5 = f.create_group('PartType5') 

part5.create_dataset('Coordinates',data=pos5)
part5.create_dataset('Velocities',data=vel5)
part5.create_dataset('ParticleIDs',data=ID5)
part5.create_dataset('Masses',data=partmass5)

f.close()


