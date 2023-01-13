import numpy as np
from scipy import constants as const

#### BULK RUN OPTION ####

run_type = 2 # 0 = save only / 1 = show only / 2 = save and show

#### RUN SETTINGS ####

N_DM = 1e7 # 1e4 / 1e5 / 1e6 / 1e7 / 1e8
M_BH = 1e-2 # 1e6 Msun / 1e7 Msun / 1e8 Msun : convert units to 10^10 Msun
R_init = np.asarray([[7, 0, 0]], dtype='float32') #BH initial position
snap_i = 100 #snap of interest 
N_PotPart = N_DM # >= N_DM

#### FIXED QUANTITIES ####

M200 = 1.0 # in units of 10^10 Msun (virial total halo mass)
R200 = 35.0 # kpc (virial radius)
V200 = 35.0 # km/s (virial velocity)
scale_length = 6.0 # kpc
h = 1 # (normalised in run)
G = const.G / 1000**3 / const.parsec * 1.989e40 # Newton's constant of gravitation in units of km^2 kpc {1e10 M_sun}^-1 s^{-2}
Gyr_sec = 1e9 * 365 * 24 * 60**2

if N_DM == 1e4:
	bin_width = 0.3
	soft_length = 0.95
elif N_DM == 1e5:
	bin_width = 0.1
	soft_length = 0.44
elif N_DM == 1e6:
	bin_width = 0.025
	soft_length = 0.20
elif N_DM == 1e7:
	bin_width = 0.01
	soft_length = 0.09
elif N_DM == 1e8:
	bin_width = 0.003
	soft_length = 0.042

b_90 = M_BH / M200 * R200 # in units of kpc
soft_length_BH = 2 * soft_length
r_sc = 2 * G * M_BH / (const.c / 1000)**2

