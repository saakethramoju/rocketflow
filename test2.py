from System import *
from Solvers import *

from constants import *


ModelNetwork = Network("Model Network")

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

SourceFluid = FluidLookup(
    "Source Fluid",
    ModelNetwork,
    {"gn2": 0.75, "O2": 0.01, "Ar": 0.24},
    pressure=3e5,
    temperature=300,
)

DarcyOption = DischargeCoefficient.model(
    "darcy",
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=101325,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A / 4,
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


Outlet1.build("darcy")

try:
    raise RuntimeError("fake darcy failure")

except Exception as e:
    print(f"{Outlet1.active_option} failed: {e}")

    Outlet1.build_next()

    SteadyState(ModelNetwork).solve(
        verbose=True,
        print_solution=True,
    )