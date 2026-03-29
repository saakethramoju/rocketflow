import numpy as np
from scipy.optimize import root_scalar

from Global import UStoSI
from Components import Source, Drain, Line, Volume

reservoir = Source(pressure=3e5)
ambient = Drain(pressure=14.67 * UStoSI.PSI_TO_PA)
manifold = Volume(pressure=1e5)
inlet = Line(cross_sectional_area=0.5e-4)
outlet = Line(cross_sectional_area=0.5e-4)

Cd = 1
rho = 1e3

def residual(P_manifold):
    dP_in = reservoir.pressure - P_manifold
    dP_out = P_manifold - ambient.pressure

    mdot_in = Cd * inlet.cross_sectional_area * np.sqrt(2 * rho * dP_in)
    mdot_out = Cd * outlet.cross_sectional_area * np.sqrt(2 * rho * dP_out)

    return mdot_in - mdot_out

# valid physical bounds: ambient < P_manifold < reservoir
P_lo = ambient.pressure + 1.0
P_hi = reservoir.pressure - 1.0

sol = root_scalar(residual, bracket=[P_lo, P_hi], method="brentq")

manifold.pressure = sol.root

print(f"Manifold Pressure: {manifold.pressure/UStoSI.PSI_TO_PA:.2f} psia")