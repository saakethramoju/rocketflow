from System import *
from Solvers import *

from constants import *


TubeNetwork = Network("Tube Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    TubeNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=304.60556,
)

ManifoldGas = IdealGasLookup(
    "Manifold gas",
    TubeNetwork,
    "gn2",
    pressure=40 * PSIA_TO_PA,
    temperature=350,
    flash_values=("pressure", "enthalpy")
)


ExitGas = IdealGasLookup(
    "Exit Gas",
    TubeNetwork,
    "gn2",
    pressure=35.78 * PSIA_TO_PA,
    temperature=417.92778,
)

L = 3207 * IN_TO_M
D = 6 * IN_TO_M
A = np.pi / 4 * D**2


Tube1 = CompressibleFlowTube(
    "Tube 1",
    TubeNetwork,
    mass_flow=10,
    upstream_static_pressure=SourceGas.pressure,
    upstream_static_temperature=SourceGas.temperature,
    upstream_density=SourceGas.density,
    downstream_static_pressure=ManifoldGas.pressure,
    downstream_static_temperature=ManifoldGas.temperature,
    downstream_density=ManifoldGas.density,
    length=L/2,
    inner_diameter=D,
    friction_factor=0,
    upstream_static_enthalpy=SourceGas.enthalpy,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    downstream_speed_of_sound=ManifoldGas.speed_of_sound,
    specific_heat_ratio=SourceGas.specific_heat_ratio
)
'''
Tube1Friction = Churchill(
    "Tube 1 Friction",
    TubeNetwork,
    mass_flow=Tube1.mass_flow,
    friction_factor=Tube1.friction_factor,
    hydraulic_diameter=D,
    dynamic_viscosity=SourceGas.dynamic_viscosity,
    cross_sectional_area=A,
    roughness=1e-6
)
'''
Manifold = Volume(
    "Manifold",
    TubeNetwork,
    pressure=ManifoldGas.pressure,
    enthalpy=ManifoldGas.enthalpy,
    volume=1,
    total_enthalpy_in=Tube1.total_enthalpy,
    heat_rate=180 * BTU_S_TO_W,
    mass_flow_in=Tube1.mass_flow,
    mass_flow_out=10
)



Tube2 = CompressibleFlowTube(
    "Tube 2",
    TubeNetwork,
    mass_flow=Manifold.mass_flow_out,
    upstream_static_pressure=ManifoldGas.pressure,
    upstream_static_temperature=ManifoldGas.temperature,
    upstream_density=ManifoldGas.density,
    downstream_static_pressure=ExitGas.pressure,
    downstream_static_temperature=ExitGas.temperature,
    downstream_density=ExitGas.density,
    length=L/2,
    inner_diameter=D,
    friction_factor=0,
    upstream_static_enthalpy=ManifoldGas.enthalpy,
    upstream_speed_of_sound=ManifoldGas.speed_of_sound,
    downstream_speed_of_sound=ExitGas.speed_of_sound,
    specific_heat_ratio=ManifoldGas.specific_heat_ratio,
    total_enthalpy=Manifold.total_enthalpy_out
)
'''

Tube2Friction = Churchill(
    "Tube 1 Friction",
    TubeNetwork,
    mass_flow=Tube2.mass_flow,
    friction_factor=Tube2.friction_factor,
    hydraulic_diameter=D,
    dynamic_viscosity=ManifoldGas.dynamic_viscosity,
    cross_sectional_area=A,
    roughness=1e-6
)
'''
solution = SteadyState(TubeNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))