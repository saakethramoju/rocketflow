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
'''
ManifoldGas = IdealGasLookup(
    "Manifold gas",
    FFNetwork,
    "gn2",
    pressure=23.4 * PSIA_TO_PA,
    temperature=288.706,
    flash_values=("pressure", "enthalpy")
)
'''

L = 3207 * IN_TO_M
D = 6 * IN_TO_M
A = np.pi / 4 * D**2




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

TubeFriction = Churchill(
    "Tube Friction",
    FFNetwork,
    mass_flow=Tube.mass_flow,
    friction_factor=Tube.friction_factor,
    hydraulic_diameter=D,
    dynamic_viscosity=SourceGas.dynamic_viscosity,
    cross_sectional_area=A,
    roughness=1e-6
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
'''


Nozzle = IsentropicAreaChange(
    "Nozzle",
    FFNetwork,
    upstream_mach_number=Tube.downstream_mach_number,
    upstream_static_pressure=SourceGas.pressure * Tube.static_pressure_ratio,
    upstream_static_temperature=SourceGas.temperature * Tube.static_temperature_ratio,
    specific_gas_constant=SourceGas.gas_constant,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    upstream_area=A,
    downstream_area=4.0 * A,
    exit_mach_regime="supersonic",
)

F = (Nozzle.mass_flow * Nozzle.downstream_mach_number * 
     np.sqrt(Nozzle.specific_heat_ratio*Nozzle.specific_gas_constant
             *Nozzle.static_temperature_ratio*Nozzle.upstream_static_temperature) +
    (Nozzle.downstream_static_pressure - 101325)*Nozzle.downstream_area)


ThrustBalance = Balance(
    "Thrust Balance",
    FFNetwork,
    variable=Nozzle.downstream_area,
    function=F - 50*LBF_TO_N
)


solution = SteadyState(FFNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))


print(f"Thrust: {F.value * N_TO_LBF: .3f} lbf")