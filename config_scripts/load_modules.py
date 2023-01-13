import numpy as np
from scipy import stats, constants as const, optimize, integrate
from scipy.special import erf
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as Rot
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import hernquist_properties as hq
from read_snap_files import *
from DF_scripts import *
from halo_calibrate import *
from config import *
import csv, time, ast
#from copy import copy

plt.rcParams["font.family"] = "Serif"
plt.rcParams["mathtext.fontset"] = "stix"
