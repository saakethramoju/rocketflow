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
D_pipe = 6 * IN_TO_M
area = np.pi / 4 * D_pipe**2




Tube = CompressibleFlowTube(
    "Tube",
    FFNetwork,
    mass_flow=1, 
    upstream_static_pressure=SourceGas.pressure,
    upstream_density=SourceGas.density,
    downstream_static_pressure=ManifoldGas.pressure,
    downstream_density=ManifoldGas.density,
    friction_factor=0.002,
    length=L,
    inner_diameter=D_pipe,

    upstream_static_enthalpy=SourceGas.enthalpy,
    upstream_speed_of_sound=SourceGas.speed_of_sound,
    downstream_speed_of_sound=ManifoldGas.speed_of_sound,

    choked=True
)

'''
TubeFriction = Churchill(
    "Tube Fanno Friction",
    FFNetwork,
    mass_flow=Tube.mass_flow,
    friction_factor=Tube.friction_factor,
    hydraulic_diameter=D_pipe,
    dynamic_viscosity=SourceGas.dynamic_viscosity,
    cross_sectional_area=area,
    roughness=1e-4,
)
'''
'''

Manifold = Volume("Intermediate",
                        FFNetwork,
                        pressure=ManifoldGas.pressure,
                        enthalpy=ManifoldGas.enthalpy,
                        volume=0.01,
                        total_enthalpy_in=Tube.total_enthalpy,
                        mass_flow_in=Tube.mass_flow)


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

mdot = Nozzle.mass_flow
Me = Nozzle.downstream_mach_number
R = Nozzle.gas_constant
k = Nozzle.specific_heat_ratio
Ti = Nozzle.upstream_static_temperature
Te_Ti = Nozzle.static_temperature_ratio
Te = Ti * Te_Ti
Pe = Nozzle.downstream_static_pressure
Ae = Nozzle.downstream_area
Pamb = 101325

F = mdot*Me*np.sqrt(k*R*Te) + (Pe - Pamb)*Ae



PeBalance = Balance("Exit Pressure Balance",
                    FFNetwork,
                    variable=Tube.length,
                    function=F - 60*LBF_TO_N)
                    #bounds=(0, None))
'''

solution = SteadyState(FFNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
)

print(solution.to_string(index=False))

#print(F.value * N_TO_LBF)
print(Tube.length.value * M_TO_IN)