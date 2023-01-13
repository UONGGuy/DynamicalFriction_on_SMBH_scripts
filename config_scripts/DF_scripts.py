import numpy as np
from scipy import constants as const
import hernquist_properties as hq
from scipy.special import erf
from halo_calibrate import *
from config import *

def X(r):
        return hq.v_circ(r) / (np.sqrt(2) * np.sqrt(3) * (hq.v_disp_radial(r)))

def decay_est_isotherm(r): # in GYr, input in kpc
        return 1.65 / np.log(np.linalg.norm(R_init)/ b_90) * np.linalg.norm(R_init)**2 * np.sqrt(3) * hq.v_disp_radial(np.linalg.norm(R_init)) / (hq.G * M_BH / const.parsec) / (3.16e16)

def decay_num_int1(r, v, rho): # fixed b_max = R_orbit, input limits in kpc, return answer in GYr
        return -(np.sqrt(G * M_BH * r) * (r + 3 * scale_length) * v**2) / (8 * np.pi * G**2 * M_BH * rho * np.log(np.linalg.norm(R_init) / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

def decay_num_int2(r, v, rho): # b_max = R_orbit(t), input limits in kpc, return answer in GYr
        return -(np.sqrt(G * M_BH * r) * (r + 3 * scale_length) * v**2) / (8 * np.pi * G**2 * M_BH * rho * np.log(r / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

def decay_num_int3(r): # b_max = rho / grad(rho), input limits in kpc, return answer in GYr, TO BE ADJUSTED
        return -(np.sqrt(G * M_BH * r) * (r + 3 * scale_length) * v_circ(r)**2) / (8 * np.pi * G**2 * M_BH * density(r) * np.log(r / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

def X_a(v, sigma):
        return np.linalg.norm(v) / (np.sqrt(2) * sigma)

def accel_est1(r, v, sigma, rho): #in units of kpc s^{-2}
        return 4 * np.pi * G**2 / v**2 * M_BH * rho * np.log(np.linalg.norm(R_init) / b_min) * (erf(X_a(v, sigma)) - 2 * X_a(v, sigma) / np.sqrt(np.pi) * np.exp(-X_a(v, sigma)**2)) / const.parsec

def accel_est_A(r, v, sigma, rho, beta): #in units of kpc s^{-2}
        return 2 * np.pi * G**2 / v**2 * M_BH * rho * np.log(1 + (np.linalg.norm(R_init) / b_min(beta))**2) * (erf(X_a(v, sigma)) - 2 * X_a(v, sigma) / np.sqrt(np.pi) * np.exp(-X_a(v, sigma)**2)) / const.parsec

def accel_est_B(r, v, sigma, rho, beta): #in units of kpc s^{-2}
	return 2 * np.pi * G**2 / np.linalg.norm(v)**3 * M_BH * rho * np.log(1 + (np.linalg.norm(R_init) / b_min(beta))**2) * (erf(X_a(v, sigma)) - 2 * X_a(v, sigma) / np.sqrt(np.pi) * np.exp(-X_a(v, sigma)**2)) / const.parsec * v


def radial_change_est(R_BH, V_BH, sigma, rho, beta): #in units of kpc s^{-1}, mutliply by dt
	return - (4 * np.pi * G**2 * rho * (R_BH + scale_length) * R_BH * M_BH * np.log(1 + (np.linalg.norm(R_init) / b_min(beta))**2)) / (V_BH**3 * (R_BH + 3 * scale_length)) / const.parsec
#	return - (4 * np.pi * G**2 * rho * (R_BH + scale_length) * R_BH * np.log(1 + (np.linalg.norm(R_init) / b_min(beta))**2)) / (V_BH**3 * (R_BH + 3 * scale_length)) * (erf(X_a(V_BH, sigma)) - 2 * X_a(V_BH, sigma) / np.sqrt(np.pi) * np.exp(-X_a(V_BH, sigma)**2)) / const.parsec


def velocity_change_est(R_BH, V_BH, sigma, rho, beta): #in units of kpc s^{-2}, multiply by dt
	return - (2 * np.pi * G**2 * rho * np.log(1 + (np.linalg.norm(R_init) / b_min(beta))**2) * M_BH) / V_BH**2 / const.parsec
#	return - (2 * np.pi * G**2 * rho * np.log(1 + (np.linalg.norm(R_init) / b_min(beta))**2)) / V_BH**2 * (erf(X_a(V_BH, sigma)) - 2 * X_a(V_BH, sigma) / np.sqrt(np.pi) * np.exp(-X_a(V_BH, sigma)**2)) / const.parsec
