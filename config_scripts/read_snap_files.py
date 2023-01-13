import numpy as np
import h5py


def get_snap_filename(input_dir, snap_no, snap_prefix="snap",
                      multiple=False, subsnap_no = 0):
    """Returns formatted file path to snap file for the respective input directory and snap number."""
    if multiple:
        snap_dir = "snapdir_%03d" %snap_no
        snap_filename = "%s/%s/%s_%03d.%d.hdf5" %(input_dir,snap_dir,
                                                  snap_prefix, snap_no,
                                                  subsnap_no)
    else:
        snap_filename = "%s/%s_%03d.hdf5" %(input_dir,snap_prefix,snap_no)
    return snap_filename


def get_attribute(fileName, attr):
    """Returns requested attribute from header of hdf5 file."""
    try:
        f = h5py.File(fileName,  "r")
    except:
        print("File %s not found" %fileName)
    h5_attr = f['/Header'].attrs[attr]
    return h5_attr


def get_snap_data(fileName, particle_type, Quant, Slicing=False, SliceInds=None):
    f = h5py.File(fileName,  "r")
    if Slicing:
        if np.asarray(SliceInds).all() == None:
            print("Need to provide slicing indices.")
        else:
            quants = f['/PartType%d/%s'%(particle_type,Quant)]
            quants = np.asarray(quants)[SliceInds]
    else:
        quants = f['/PartType%d/%s'%(particle_type,Quant)]
    return np.asarray(quants)


def get_snap_attr_2(fileName, snap_no, attr):
	f = h5py.File(fileName, 'r')
	h5_attr = f['/snap_%03d' %(snap_no)].attrs[attr]
	return h5_attr


def get_snap_data_2(fileName, snap_no, Quant):
	f = h5py.File(fileName, 'r')
	quants = f['/snap_%03d/%s/' %(snap_no, Quant)]
	return np.asarray(quants)


