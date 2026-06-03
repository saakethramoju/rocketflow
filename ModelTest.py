from System import *
from Solvers import *

from constants import *


ModelNetwork = Network("Multiple Model Test")

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

SourceFluid = FluidLookup(
    "Source Fluid",
    ModelNetwork,
    "rp-1",
    pressure=3e5,
    temperature=300,
)

# ------------------------------------------------------------------
# Model 1: Outlet 1
# First option is intentionally bad.
# ------------------------------------------------------------------

Outlet1MassFlow = State(1.0)
Outlet1FrictionFactor = State(0.02)
BadDensity = State(-1.0)

Outlet1BadDarcy = ModelOption(
    "bad_darcy",
    components=[
        Colebrook.model(
            "Outlet 1 Bad Friction",
            mass_flow=Outlet1MassFlow,
            friction_factor=Outlet1FrictionFactor,
            hydraulic_diameter=D,
            dynamic_viscosity=SourceFluid.dynamic_viscosity,
            cross_sectional_area=A / 4,
        ),
        DarcyWeisbach.model(
            "Outlet 1 Bad Darcy",
            mass_flow=Outlet1MassFlow,
            upstream_pressure=SourceFluid.pressure,
            downstream_pressure=101325,
            length=1,
            cross_sectional_area=A / 4,
            hydraulic_diameter=D,
            density=BadDensity,
            friction_factor=Outlet1FrictionFactor,
        ),
    ],
)

Outlet1CdA = DischargeCoefficient.model(
    "cda",
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=101325,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A / 4,
)

Outlet1 = Model(
    "Outlet 1",
    ModelNetwork,
    Outlet1BadDarcy,
    Outlet1CdA,
    order=[
        "bad_darcy",
        "cda",
    ],
)

# ------------------------------------------------------------------
# Model 2: Outlet 2
# Multiple valid options.
# Only the first/default option should be used during Outlet 1 sweep.
# ------------------------------------------------------------------

Outlet2MassFlow = State(1.0)
Outlet2FrictionFactor = State(0.02)

Outlet2Darcy = ModelOption(
    "darcy",
    components=[
        Churchill.model(
            "Outlet 2 Friction",
            mass_flow=Outlet2MassFlow,
            friction_factor=Outlet2FrictionFactor,
            hydraulic_diameter=D,
            dynamic_viscosity=SourceFluid.dynamic_viscosity,
            cross_sectional_area=A,
        ),
        DarcyWeisbach.model(
            "Outlet 2 Darcy",
            mass_flow=Outlet2MassFlow,
            upstream_pressure=SourceFluid.pressure,
            downstream_pressure=101325,
            length=1,
            cross_sectional_area=A,
            hydraulic_diameter=D,
            density=SourceFluid.density,
            friction_factor=Outlet2FrictionFactor,
        ),
    ],
)

Outlet2CdA = DischargeCoefficient.model(
    "cda",
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=101325,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
)

Outlet2 = Model(
    "Outlet 2",
    ModelNetwork,
    Outlet2Darcy,
    Outlet2CdA,
    order=[
        "darcy",
        "cda",
    ],
)

# ------------------------------------------------------------------
# Sweep Outlet 1.
# Expected:
#   - Outlet 1 bad_darcy fails and is skipped.
#   - Outlet 1 cda succeeds.
#   - Outlet 2 has multiple options, but only its default darcy option is used.
# ------------------------------------------------------------------

results = SteadyState(ModelNetwork).solve(
    model="Outlet 1",
    evaluate_all_model_options=True,
    verbose=True,
    print_solution=True,
    filename="test.xlsx",
)