import numpy as np
from scipy.optimize import root_scalar

from Components import Source, Drain

s = Drain("tank", mass_flow=3)
print(s)



