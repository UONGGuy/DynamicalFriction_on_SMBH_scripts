import numpy as np
from scipy import stats
from config import *

def b_min(beta):
        return np.sqrt(b_90**2 + beta * soft_length_BH**2 + r_sc**2)

def find_CoM(R, M):
        R_CoM = np.average(R, axis=0, weights=M)
        return R_CoM

def dist_CoM(CoM, R_in):
        R1 = R_in - CoM
        Mag = np.linalg.norm(R1, axis=1)
        R2 = np.c_[R1, Mag]
        return R2

def radial_density_halo(D_CoM, M, dr, R200):
        M_binned, bin_edges, no_bins = stats.binned_statistic(D_CoM, M, 'sum', np.arange(0, R200 * 1.1, dr))
        dVol = 4 * np.pi * (bin_edges[1:] ** 2) * dr
        rho = np.divide(M_binned, dVol)
        R = np.column_stack((bin_edges[1:], rho))
        return R
