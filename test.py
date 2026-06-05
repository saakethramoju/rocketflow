from System import *
from Solvers import *
from constants import *

ThermalSystem = Network("Thermal")

D = 0.5 * IN_TO_M
L = 10.0 * IN_TO_M

flow_area = (np.pi / 4.0) * D**2
wetted_area = np.pi * D * L

TubeMaterial = MaterialLookup(
    "Tube Material",
    ThermalSystem,
    "718",
)

SourceFluid = FluidLookup(
    "Fluid Source",
    ThermalSystem,
    "lox",
    pressure=5e5,
    quality=0,
)

VolFluid = FluidLookup(
    "Fluid Volume",
    ThermalSystem,
    "lox",
    pressure=1e5,
    quality=0,
)

TubeFlow = DarcyWeisbach(
    "Tube Flow",
    ThermalSystem,
    mass_flow=5,
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=VolFluid.pressure,
    length=L,
    cross_sectional_area=flow_area,
    hydraulic_diameter=D,
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
    roughness=1e-9,
)

Tube = Solid(
    "Tube",
    ThermalSystem,
    temperature=TubeMaterial.temperature,
)

LOXConvection = Convection(
    "LOX Convection",
    ThermalSystem,
    surface_temperature=Tube.temperature,
    fluid_temperature=SourceFluid.temperature,
    convective_area=wetted_area,
    convection_coefficient=25,
)

LOXDittus = DittusBoelter(
    "LOX Dittus",
    ThermalSystem,
    hydraulic_diameter=TubeFlow.hydraulic_diameter,
    fluid_conductivity=SourceFluid.conductivity,
    fluid_specific_heat=SourceFluid.specific_heat_cp,
    fluid_dynamic_viscosity=SourceFluid.dynamic_viscosity,
    cross_sectional_area=TubeFlow.cross_sectional_area,
    mass_flow=TubeFlow.mass_flow,
    convection_coefficient=LOXConvection.convection_coefficient,
)

AirConvection = Convection(
    "Air Convection",
    ThermalSystem,
    surface_temperature=Tube.temperature,
    fluid_temperature=293.15,
    convective_area=wetted_area,
    convection_coefficient=25,
)

Tube.heat_rate = LOXConvection.heat_rate + AirConvection.heat_rate

solution = SteadyState(ThermalSystem).solve(
    verbose=True,
    print_solution=True,
)