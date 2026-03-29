import numpy as np

from Components import Source, Drain, Line

s = Source(pressure=2e5)

d = Drain(pressure=101325)

b = Line(cross_sectional_area=0.5e-4, length=1)

Cd = 1
rho = 1e3

dp = s.pressure - d.pressure
b.mass_flow = Cd * b.cross_sectional_area * np.sign(dp) * np.sqrt(2 * rho * abs(dp))
print(b)