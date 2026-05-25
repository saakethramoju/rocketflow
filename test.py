from System import *
from Solvers import *

from constants import *


FFNetwork = Network("Fanno Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    FFNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=299.817,
)

ManifoldGas = IdealGasLookup(
    "Manifold gas",
    FFNetwork,
    "gn2",
    pressure=23.4 * PSIA_TO_PA,
    temperature=288.706,
    flash_values=("pressure", "enthalpy")
)


L = 3207 * IN_TO_M
D = 6 * IN_TO_M
area = np.pi / 4 * D**2



'''
Tube = ChokedFannoFlow(
    "Tube",
    FFNetwork,
    upstream_density=SourceGas.density,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    friction_factor=0.002,
    length=L,
    inner_diameter=D,
    upstream_static_enthalpy=SourceGas.enthalpy,
    regime="subsonic"
)
'''

Tube = ChokedFannoFlow(
    "Tube",
    FFNetwork,
    upstream_density=SourceGas.density,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    friction_factor=0.002,
    length=L,
    inner_diameter=D,
    upstream_static_enthalpy=SourceGas.enthalpy,
    upstream_mach_number=1.3,
    regime="supersonic",
)

solution = SteadyState(FFNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))

