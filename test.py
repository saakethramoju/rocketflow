from System import *
from Solvers import *
from constants import *
from matprotlib import Steel

T = 300
steel = Steel("1018 Carbon")
k = steel.get("thermal_conductivity", T)
cp = steel.get("specific_heat", T)

L = 0.5 * IN_TO_M
D = 1 * IN_TO_M
A = (np.pi/4) * D**2

ThermalSystem = Network("Thermal")


NodeTemp = State(300)

ThermalSystem.track(
    "Node Temp",
    NodeTemp
)

Rod1 = SolidConductor(
    "Rod 1",
    ThermalSystem,
    upstream_temperature=500,
    downstream_temperature=NodeTemp,
    thermal_conductivity=k,
    length=L,
    cross_sectional_area=A,
)

Metal = Solid(
    "Solid Node",
    ThermalSystem,
    temperature=Rod1.downstream_temperature,
    heat_rate_in=Rod1.heat_rate
)

Rod2 = SolidConductor(
    "Rod 2",
    ThermalSystem,
    upstream_temperature=Metal.temperature,
    downstream_temperature=300,
    thermal_conductivity=k,
    length=L,
    cross_sectional_area=A,
    heat_rate=Metal.heat_rate_out
)

solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)
