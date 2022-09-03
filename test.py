from ctypes import * 
from numpy.ctypeslib import ndpointer

lib = CDLL(r"test.dll", winmode=0)
lib.simulate.restype = ndpointer(dtype=c_double, shape=(5000,))

zeroes = (c_double * 1)(-1/5,)
poles = (c_double * 2)(0, -1/3)
overall_gain = c_double(-250 * 65536 / 6)

# check if size matches 
if lib.get_len_zeroes() != len(zeroes):
    raise ValueError('Size of zeroes differ')
if lib.get_len_poles() != len(poles):
    raise ValueError('Size of poles differ')

from time import time 
tt = time() 
for _ in range(10000):
    errs = lib.simulate(zeroes, poles, overall_gain)
print(time() - tt)
# compare the results with Python 
from cascade_sim import IIRCascadeRealSimulator
from sim_free import noise
def plant(dac_num):
    return 4.1e-4  * dac_num * 6 / 65536

zeroes, poles, gain  = [-1/5], [0, -1/3], -250 # specs for analog filter 
iir = IIRCascadeRealSimulator(zeroes, poles, gain * 65536 / 6, plant)  # define 
errors = []
for e in noise[:5000]: 
    errors.append(iir.update(e))

# should print False, otherwise it's an error 
print(any(e_cpp != e_py for e_cpp, e_py in zip(errs, errors)))

