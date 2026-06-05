from System import *
from Solvers import *
from constants import *


ThermalSystem = Network("Thermal")

mat = "c17200"
L = 0.5*IN_TO_M
D = 1 * IN_TO_M
A = (np.pi/4) * D**2


SourceMaterial = MaterialLookup(
    "Source Material",
    ThermalSystem,
    mat,
    temperature=600
)


NodeMaterial = MaterialLookup(
    "Node Material",
    ThermalSystem,
    mat,
)



Rod1 = SolidConductor(
    "Rod 1",
    ThermalSystem,
    upstream_temperature=SourceMaterial.temperature,
    downstream_temperature=NodeMaterial.temperature,
    thermal_conductivity=SourceMaterial.thermal_conductivity,
    length=L,
    cross_sectional_area=A,
)


Metal = Solid(
    "Solid Node",
    ThermalSystem,
    temperature=NodeMaterial.temperature,
    heat_rate_in=Rod1.heat_rate
)


Rod2 = SolidConductor(
    "Rod 2",
    ThermalSystem,
    upstream_temperature=Metal.temperature,
    downstream_temperature=300,
    thermal_conductivity=NodeMaterial.thermal_conductivity,
    length=L,
    cross_sectional_area=A,
    heat_rate=Metal.heat_rate_out
)


solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)
