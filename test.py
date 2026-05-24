from System import *
from Solvers import *

from constants import *


FFNetwork = Network("Fanno Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    FFNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=3000,
)

ManifoldGas = IdealGasLookup(
    "Manifold gas",
    FFNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=3000,
    flash_values=("pressure", "enthalpy")
)


L = 6 * IN_TO_M
D = 2.875 * IN_TO_M
area = np.pi / 4 * D**2

'''
Tube = ChokedFannoFlow(
    "Fanno Tube",
    FFNetwork,
    upstream_density=SourceGas.density,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    friction_factor=0.002,
    length=L,
    inner_diameter=D,
    regime="subsonic",
    upstream_static_enthalpy=SourceGas.enthalpy,

)
'''


Tube = UnchokedFannoFlow(
    "Tube",
    FFNetwork,
    upstream_density=SourceGas.density,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    downstream_density=ManifoldGas.density,
    downstream_speed_of_sound=ManifoldGas.speed_of_sound,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    friction_factor=0.002,
    length=L,
    inner_diameter=D,
    upstream_static_enthalpy=SourceGas.enthalpy
)


TubeFriction = Churchill(
    "Tube Fanno Friction",
    FFNetwork,
    mass_flow=Tube.mass_flow,
    friction_factor=Tube.friction_factor,
    hydraulic_diameter=D,
    dynamic_viscosity=SourceGas.dynamic_viscosity,
    cross_sectional_area=area,
    roughness=1e-4,
)



Manifold = Volume("Intermediate",
                        FFNetwork,
                        pressure=ManifoldGas.pressure,
                        enthalpy=ManifoldGas.enthalpy,
                        volume=0.01,
                        total_enthalpy_in=Tube.total_enthalpy,
                        mass_flow_in=Tube.mass_flow)

'''
Manifold = IsothermalVolume("Intermediate",
                            FFNetwork,
                            pressure=ManifoldGas.pressure,
                            temperature=ManifoldGas.temperature,
                            volume=0.01,
                            mass_flow_in=Tube.mass_flow)
'''

Nozzle = IsentropicAreaChange(
    "Nozzle",
    FFNetwork,
    upstream_mach_number=Tube.downstream_mach_number,
    upstream_static_pressure=Manifold.pressure,
    upstream_static_temperature=ManifoldGas.temperature,
    gas_constant=ManifoldGas.gas_constant,
    specific_heat_ratio=ManifoldGas.specific_heat_ratio,
    upstream_area=area,
    downstream_area=1.1*area,
    mass_flow=Manifold.mass_flow_out,
    total_enthalpy=Manifold.total_enthalpy_out,
    exit_mach_regime='supersonic'
)

'''
PeBalance = Balance("Exit Pressure Balance",
                    FFNetwork,
                    variable=Tube.length,
                    function=Nozzle.downstream_static_pressure - 2*PSIA_TO_PA)
'''

solution = SteadyState(FFNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))