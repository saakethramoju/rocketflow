from System import *
from Solvers import *

from constants import *


TubeNetwork = Network("Tube Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    TubeNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=299.817,
)

ManifoldGas = IdealGasLookup(
    "Manifold gas",
    TubeNetwork,
    "gn2",
    pressure=40 * PSIA_TO_PA,
    temperature=350,
    flash_values=("pressure", "enthalpy")
)


L = 3207 * IN_TO_M
D = 0.25 * IN_TO_M
A = np.pi / 4 * D**2

Tube = ChokedRayleighFlow(
    "Tube",
    TubeNetwork,
    upstream_density=SourceGas.density,
    upstream_static_temperature=SourceGas.temperature,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    specific_gas_constant=SourceGas.gas_constant,
    inner_diameter=D,
    heat_rate=1800 * BTU_S_TO_W,
    upstream_static_enthalpy=SourceGas.enthalpy,
    regime="supersonic",
    #upstream_mach_number=4
)

solution = SteadyState(TubeNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))