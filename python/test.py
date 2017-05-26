
import tetrachotomy
import importlib
importlib.reload(tetrachotomy)


import time
import tetrachotomy as tc


tetra = tc.PoleHunt()

def func(z):
    time.sleep(1e-1)
    return 1/z

integral, I = tetra.compute_integral(func, 0)
