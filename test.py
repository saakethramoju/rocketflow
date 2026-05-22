from System import *
from Solvers import *

from constants import *

FFNetwork = Network("Fanno Flow")

Source = IdealGasLookup(
            "Fuel Tank Ullage Gas",
            FFNetwork,
            'gn2',
            pressure=50 * PSIA_TO_PA,
            temperature=300,
        )

Drain = IdealGasLookup(
    "atmosphere",
    FFNetwork,
    'gn2',
    pressure=23.4*PSIA_TO_PA,
    temperature=300
)

a = Source.speed_of_sound
rho = Source.density
D = 6 *IN_TO_M
area = np.pi/4 * (D)**2

mdot_0 = 0.5*rho*a*area


Tube = FannoFlow("Fanno Tube",
                 FFNetwork,
                 mass_flow=mdot_0.value,
                 upstream_pressure=Source.pressure,
                 upstream_density=Source.density,
                 upstream_speed_of_sound=Source.speed_of_sound,
                 upstream_specific_heat_ratio=Source.specific_heat_ratio,
                 downstream_pressure=Drain.pressure,
                 downstream_density=Drain.density,
                 downstream_speed_of_sound=Drain.speed_of_sound,
                 downstream_specific_heat_ratio=Drain.specific_heat_ratio,
                 length=3207 * IN_TO_M,
                 inner_diameter=6 * IN_TO_M,
                 friction_factor=0.02)


TubeFriction = Churchill("Tube Fanno Friction",
                         FFNetwork,
                         mass_flow=Tube.mass_flow,
                         friction_factor=Tube.friction_factor,
                         hydraulic_diameter=D,
                         dynamic_viscosity=Source.dynamic_viscosity,
                         cross_sectional_area=area,
                         roughness=0.1e-4)

solution = SteadyState(FFNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))

