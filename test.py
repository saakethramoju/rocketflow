from System import *
from Solvers import *

from constants import *


FFNetwork = Network("Fanno Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    FFNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=300,
)

ManifoldGas = IdealGasLookup(
    "Manifold gas",
    FFNetwork,
    "gn2",
    pressure=23.4 * PSIA_TO_PA,
    temperature=288.706,
    flash_values=("pressure", "enthalpy")
)

L = 10 * IN_TO_M
D = 6 * IN_TO_M
area = np.pi / 4 * D**2


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

k = ManifoldGas.specific_heat_ratio
Po = ManifoldGas.pressure * (1 + (k-1)/2)**(k/(k-1))
To = ManifoldGas.temperature * (1 + (k-1)/2)

Orifice = IsentropicCompressibleOrifice("Orifice",
                                        FFNetwork,
                                        upstream_total_pressure=Po,
                                        upstream_total_temperature=To,
                                        downstream_pressure=101325,
                                        discharge_coefficient=1,
                                        cross_sectional_area=(np.pi/4)*(1*IN_TO_M)**2,
                                        specific_gas_constant=ManifoldGas.gas_constant,
                                        specific_heat_ratio=ManifoldGas.specific_heat_ratio,
                                        mass_flow=Manifold.mass_flow_out,
                                        total_enthalpy=Manifold.total_enthalpy_out)




solution = SteadyState(FFNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))
