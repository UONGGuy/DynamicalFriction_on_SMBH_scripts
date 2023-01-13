from load_modules import *
from config import *

print(b_90, r_sc)
print(b_min(0), b_min(1), b_min(1.5))
print(np.log(1 + (7/b_min(0))**2), np.log(1 + (7/b_min(1))**2), np.log(1 + (7/b_min(1.5))**2))
