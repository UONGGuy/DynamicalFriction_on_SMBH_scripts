from load_modules import *

def X(r):
	return hq.v_circ(r) / (np.sqrt(2) * np.sqrt(3) * (hq.v_disp_radial(r)))

def r_dot(r):
        return ((- 8 * np.pi * hq.G**2 * r * (r + scale_length)**2) / (np.sqrt(hq.G * M_BH * r) * (3 * scale_length + r) * hq.v_circ(r)**2) * M_BH * hq.density(r) * np.log(np.linalg.norm(R_init) / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2))) / const.parsec

def isotherm(r):
	return r / (-0.302 *  np.log(np.linalg.norm(R_init) / b_90) * hq.G / const.parsec * M_BH / (np.sqrt(3) * hq.v_disp_radial(np.linalg.norm(R_init))))

def integrand(r):
	return -(np.sqrt(hq.G * M_BH* r) * (r + 3 * scale_length) * hq.v_circ(r)**2) / (8 * np.pi * hq.G**2 * M_BH * hq.density(r) * np.log(np.linalg.norm(R_init) / b_min) * (erf(X(r)) - 2 * X(r) / np.sqrt(np.pi) * np.exp(-X(r)**2)) * r * (r + scale_length)**2) * const.parsec / (3.16e16)

#print(r_dot(R200/2 * const.parsec))
#print(r_dot(2.8 * soft_length * const.parsec))
#print((R200/2 - 2.8 * soft_length) * const.parsec / r_dot(R200 * const.parsec) / (3.16e16))

#x = lambda r : 1/r_dot(r, np.linalg.norm(R_init))
x = lambda r : 1/r_dot(r)

t, err = integrate.quad(x, np.linalg.norm(R_init), 2.8 * soft_length)

print('Analytic num int 1 (b_max = const) = ', t / (3.16e16))

t2, err2 = integrate.quad(integrand, np.linalg.norm(R_init), 2.8 * soft_length)

print('Analytic num int 2 (b_max = const) = ', t2)

ti, erri = integrate.quad(isotherm, np.linalg.norm(R_init),  2.8 * soft_length)

print('Isotherm integrated = ', ti / (3.16e16))

T1, Err1 = integrate.quad(hq.decay_est_int1, np.linalg.norm(R_init), 2.8 * soft_length)

print('HQ (const b_max) =', T1)

T2, Err2 = integrate.quad(hq.decay_est_int2, np.linalg.norm(R_init), 2.8 * soft_length)

print('HQ (b_max = R(t)) =', T2)

print('Isotherm analytic = ', decay_est_isotherm(np.linalg.norm(R_init)))

print('b_90 =', b_90)

print('b_min =', b_min)
