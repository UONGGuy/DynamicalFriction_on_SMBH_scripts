import numpy as np
from scipy import stats, constants as const
from scipy.special import erf
from config import *

#### FUNCTIONS TO RETURN VALUES FOR SPHERICAL ISOTROPIC HERNQUIST MODEL ####

# see Hernquist 1990 for more details 

a = scale_length # in units of kpc
#defined in config: G = const.G / 1000**3 / const.parsec * 1.989e40 # Newton's constant of gravitation in units of km^2 kpc {1e10 M_sun}^-1 s^{-2}
r_half_mass = (1 + np.sqrt(2)) * a # in units of kpc
kinetic_energy_tot = - G * M200**2 / (6 * a) # in units of km^2 s^{-2} 
gamma = 1 # constant mass-to-light ratio
v_g = np.sqrt((G * M200) / a) # in units of km s^1

def density(r): # in units of {1e10 M_sun} kpc^{-3}
        rho = M200 / (2.0 * np.pi) * (a / r) * 1.0 / (r + a) ** 3.0
        return rho

def cum_mass(r): # in units of {1e10 M_sun}
	M_cum = M200 * r**2 / ((r + a)**2)
	return M_cum

def potential(r): # in units of km^2 s^{-2} 
	phi = - (G * M200) / (r + a)
	return phi

def density_dimless(r):
	rho_dl = 2 * np.pi * a**3 / M200 * density(r)
	return rho_dl

def potential_dimless(r):
	phi_dl = a / (r + a)
	return phi_dl

def v_disp_radial(r): # in units of km s^{-1}, input in kpc 
	v_r2 = G * M200 / (12 * a) * ((12 * r * (r + a)**3) / a**4 * np.log((r + a) / r) - r / (r + a) * (25 + 52 * r / a + 42 * (r / a)**2 + 12 * (r / a)**3))
	return np.sqrt(v_r2)

def kinetic_energy(r): # in units of km^3 s^{-2}
	KE = const.G * M200**2 / (4 * a) * (4 * (r / a)**3 * np.log((r + a)/ r) - 4 * (r / a)**2 + 2 * r / a - 1 + (((r / a)**2 + r / a + 1) / (1 + (r / a)**3)))
	return KE

def v_esc(r): # in units of km s^{-1}
	v_e = np.sqrt((2 * G * M200) / (r + a))
	return v_e

def v_circ(r): # in units of km s^{-1}
	v_c = np.sqrt((G * M200 * r)) / (r + a)
	return v_c

def X(s):
	for element in s:
		if element < 0:
			print('\nScaled projected radius negative - should be >= 0. Element ' + element)
		if element > 1:
			return 1 / np.sqrt(s**2 - 1) * np.arccos(1 /s)
		if element <= 1:
			return 1 / np.sqrt(1 - s**2) * np.log((1 + np.sqrt(1 - s**2)) / s)

def I(R): # surface brightness function. R is projected radius. in units of {1e10 M_sun} kpc^{-2} gamma^{-1}
	s = R / a
	return M200 / (2 * np.pi * a**2 * gamma * (1 - s**2)**2) * ((2 + s**2) * X(s) - 3)

def v_disp_los(R): # if model isotropic. in units of km s^{-2}
	s = R / a
	sigma_p2 = G * M200**2 / (12 * np.pi * a**3 * I(R) * gamma) * (0.5 * (1 - s**2)**(-3) * (-3 * s**2 * X(s) * (8 * s**6 - 28 * s**4 + 35 * s**2 - 20) - 24 * s**6 + 68 * s**4 - 65 * s**2 + 6) - 6 * np.pi * s)
	return sigma_p2

def X(r): # for orbital decay est.
        return v_circ(r) / (np.sqrt(2) * np.sqrt(3) * (v_disp_radial(r)))

def decay_est_int1(r): # fixed b_max = R_orbit, input limits in kpc, return answer in GYr
        return -(np.sqrt(G * M_BH * r) * (r + 3 * scale_length) * v_circ(r)**2) / (8 * np.pi * G**2 * M_BH * density(r) * np.log(np.linalg.norm(R_init) / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

def decay_est_int2(r): # b_max = R_orbit(t), input limits in kpc, return answer in GYr
        return -(np.sqrt(G * M_BH * r) * (r + 3 * scale_length) * v_circ(r)**2) / (8 * np.pi * G**2 * M_BH * density(r) * np.log(r / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

def decay_est_int3(r): # b_max = rho / grad(rho), input limits in kpc, return answer in GYr, TO BE ADJUSTED
        return -(np.sqrt(G * M_BH * r) * (r + 3 * scale_length) * v_circ(r)**2) / (8 * np.pi * G**2 * M_BH * density(r) * np.log(r / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

