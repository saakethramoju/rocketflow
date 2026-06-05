from System import *
from Solvers import *
from constants import *


ThermalSystem = Network("Thermal")


def harmonic_mean(a, b):
    return 2*a*b/(a+b)



TubeMaterial = MaterialLookup(
    "Node Material",
    ThermalSystem,
    "718",
)

SourceFluid = FluidLookup(
    "Fluid Source",
    ThermalSystem,
    "lox",
    pressure=5e5,
    quality=0
)

VolFluid = FluidLookup(
    "Fluid Source",
    ThermalSystem,
    "lox",
    pressure=1e5,
    quality=0
)






TubeFlow = DarcyWeisbach(
    "Tub Flow",
    ThermalSystem,
    mass_flow=5,
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=VolFluid.pressure,
    length=10 * IN_TO_M,
    cross_sectional_area=(np.pi/4) * (0.5*IN_TO_M)**2,
    hydraulic_diameter=0.5*IN_TO_M,
    density=SourceFluid.density,
    friction_factor=2e-5,
)

TubeFriction = Colebrook(
    "Tube Friction",
    ThermalSystem,
    mass_flow=TubeFlow.mass_flow,
    friction_factor=TubeFlow.friction_factor,
    hydraulic_diameter=TubeFlow.hydraulic_diameter,
    dynamic_viscosity=SourceFluid.dynamic_viscosity,
    cross_sectional_area=TubeFlow.cross_sectional_area,
    roughness=1e-9
)

Tube = Solid(
    "Tube",
    ThermalSystem,
    temperature=TubeMaterial.temperature,
)

FlowConvection = Convection(
    "Flow Convection",
    ThermalSystem,
    surface_temperature=Tube.temperature,
    fluid_temperature=SourceFluid.temperature,
    convective_area=2 * TubeFlow.cross_sectional_area,
    heat_rate=Tube.heat_rate_in,
)

FlowDittus = DittusBoelter(
    "Flow Dittus",
    ThermalSystem,
    hydraulic_diameter=TubeFlow.hydraulic_diameter,
    fluid_conductivity=SourceFluid.conductivity,
    fluid_specific_heat=SourceFluid.specific_heat_cp,
    fluid_dynamic_viscosity=SourceFluid.dynamic_viscosity,
    cross_sectional_area=TubeFlow.cross_sectional_area,
    mass_flow=TubeFlow.mass_flow,
    convection_coefficient=FlowConvection.convection_coefficient
)

AirConvection = Convection(
    "Flow Convection",
    ThermalSystem,
    surface_temperature=Tube.temperature,
    fluid_temperature=293.15,
    convective_area=2 * TubeFlow.cross_sectional_area,
    convection_coefficient=25,
    heat_rate=Tube.heat_rate_out,
)

solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)
