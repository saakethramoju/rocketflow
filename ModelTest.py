from System import *
from Solvers import *

from constants import *



def solve_with_model_fallback(
    network,
    models: list,
    solver_class,
    **solve_kwargs,
):
    """
    Try solving a network while advancing model options after failures.

    This function assumes each Model has already been built or can build its
    first option. If the solve fails, the first Model with another available
    option is advanced and the solve is tried again.
    """

    # Build any unbuilt models using their first option.
    for model in models:
        if model.active_component is None:
            model.build()

    while True:
        try:
            return solver_class(network).solve(**solve_kwargs)

        except Exception as error:
            # Find the first model that still has another option.
            fallback_model = next(
                (model for model in models if model.has_next),
                None,
            )

            if fallback_model is None:
                raise RuntimeError(
                    "All model fallback options failed."
                ) from error

            print(
                f"{fallback_model.name}: "
                f"{fallback_model.active_option} failed. "
                f"Trying {fallback_model.next()}."
            )

            fallback_model.build_next()


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


solution = solve_with_model_fallback(
    ModelNetwork,
    models=[
        Outlet1,
    ],
    solver_class=SteadyState,
    verbose=True,
    print_solution=True,
)