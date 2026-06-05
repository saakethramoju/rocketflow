from System import *
from Solvers import *
from constants import *


ThermalSystem = Network("Thermal")

mat = "c17200"
L = 0.5*IN_TO_M
D = 1 * IN_TO_M
A = (np.pi/4) * D**2

def geometric_mean(a, b):
    return 2*a*b/(a+b)

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

k1 = geometric_mean(SourceMaterial.thermal_conductivity, NodeMaterial.thermal_conductivity)

Rod1 = Conduction(
    "Rod 1",
    ThermalSystem,
    temperature1=SourceMaterial.temperature,
    temperature2=NodeMaterial.temperature,
    thermal_conductivity=k1,
    length=L,
    conductive_area=A,
)


Metal = Solid(
    "Solid Node",
    ThermalSystem,
    temperature=NodeMaterial.temperature,
    heat_rate_in=Rod1.heat_rate
)


Rod2 = Conduction(
    "Rod 2",
    ThermalSystem,
    temperature1=Metal.temperature,
    temperature2=300,
    thermal_conductivity=NodeMaterial.thermal_conductivity,
    length=L,
    conductive_area=A,
    heat_rate=Metal.heat_rate_out
)


solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)
