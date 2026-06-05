from System import *
from Solvers import *
from constants import *


ThermalSystem = Network("Thermal")

mat = "718"
L = 0.5*IN_TO_M
D = 1 * IN_TO_M
A = (np.pi/4) * D**2

def harmonic_mean(a, b):
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

k1 = harmonic_mean(SourceMaterial.thermal_conductivity, NodeMaterial.thermal_conductivity)

Conductor1 = Conduction(
    "Conductor 1",
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
    heat_rate_in=Conductor1.heat_rate
)

AirConvection = Convection(
    "Air Convection",
    ThermalSystem,
    surface_temperature=Metal.temperature,
    fluid_temperature=298.15,
    heat_transfer_coefficient=20,
    convective_area=A*6
)


Conductor2 = Conduction(
    "Conductor 2",
    ThermalSystem,
    temperature1=Metal.temperature,
    temperature2=300,
    thermal_conductivity=NodeMaterial.thermal_conductivity,
    length=L,
    conductive_area=A/4,
    #heat_rate=Metal.heat_rate_out
)

Radiator1 = AmbientRadiation(
    "Radiator 1",
    ThermalSystem,
    solid_temperature=Metal.temperature,
    ambient_temperature=298.15,
    emissivity=0.8,
    radiative_area=A,
)

Metal.heat_rate_out = Conductor2.heat_rate + Radiator1.heat_rate + AirConvection.heat_rate

solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)
