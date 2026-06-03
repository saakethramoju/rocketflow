from System import *
from Solvers import *

from constants import *


ModelNetwork = Network("Model Network")

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

# ------------------------------------------------------------------
# Fluid lookup
# ------------------------------------------------------------------

SourceFluid = FluidLookup(
    "Source Fluid",
    ModelNetwork,
    "rp-1",
    pressure=3e5,
    temperature=300,
)

# ------------------------------------------------------------------
# Shared states for the Darcy/Colebrook option
# ------------------------------------------------------------------

DarcyMassFlow = State(1.0)
DarcyFrictionFactor = State(0.02)

# ------------------------------------------------------------------
# Model options
# ------------------------------------------------------------------

DarcyOption = ModelOption(
    "darcy",
    components=[
        Colebrook.model(
            "Outlet 1 Friction",
            mass_flow=DarcyMassFlow,
            friction_factor=DarcyFrictionFactor,
            hydraulic_diameter=D,
            dynamic_viscosity=SourceFluid.dynamic_viscosity,
            cross_sectional_area=A,
        ),
        DarcyWeisbach.model(
            "Outlet 1 Darcy",
            mass_flow=DarcyMassFlow,
            upstream_pressure=SourceFluid.pressure,
            downstream_pressure=101325,
            length=1,
            cross_sectional_area=A,
            hydraulic_diameter=D,
            density=SourceFluid.density,
            friction_factor=DarcyFrictionFactor,
        ),
    ],
)

CompressibleOption = CompressibleFlowTube.model(
    "compressible",
    mass_flow=State(1.0),
    upstream_static_pressure=SourceFluid.pressure,
    upstream_static_temperature=SourceFluid.temperature,
    upstream_density=SourceFluid.density,
    downstream_static_pressure=101325,
    downstream_static_temperature=SourceFluid.temperature,
    downstream_density=SourceFluid.density,
    length=1,
    inner_diameter=D,
    friction_factor=0.02,
    upstream_static_enthalpy=SourceFluid.enthalpy,
    upstream_speed_of_sound=SourceFluid.speed_of_sound,
    specific_heat_ratio=SourceFluid.specific_heat_ratio,
)

Outlet1 = Model(
    "Outlet 1",
    ModelNetwork,
    components=[
        DarcyOption,
        CompressibleOption,
    ],
    order=[
        "darcy",
        "compressible",
    ],
)


results = SteadyState(ModelNetwork).solve(
    model="Outlet 1",
    verbose=True,
    print_solution=True,
    filename="test.xlsx",
    #static=True,
    evaluate_all_model_options=True
)