from System import *
from Solvers import *

from constants import *


ModelNetwork = Network("Model Network")

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

SourceFluid = FluidLookup(
    "Source Fluid",
    ModelNetwork,
    {"gn2": 0.75, "O2": 0.01, "Ar": 0.24},
    pressure=3e5,
    temperature=300,
)

DarcyOption = DischargeCoefficient.model(
    "darcy",
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=101325,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A / 4,
)

DarcyComponent = DarcyOption.build(
    "Outlet 1",
    ModelNetwork,
)

print(DarcyComponent)
print(ModelNetwork.components)