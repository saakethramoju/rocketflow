from System import *
from Solvers import *

from constants import *


MixtureNetwork = Network("Mixture Flow")

# ------------------------------------------------------------------
# Geometry
# ------------------------------------------------------------------

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

# ------------------------------------------------------------------
# Shared outlet compositions
# ------------------------------------------------------------------

SeparatorOutlet1Composition = Composition({"Ar": 1.0})
SeparatorOutlet2Composition = Composition()

# ------------------------------------------------------------------
# Fluid lookups first
# ------------------------------------------------------------------

SourceFluid = FluidLookup(
    "Source Fluid",
    MixtureNetwork,
    {"gn2": 0.75, "O2": 0.01, "Ar": 0.24},
    pressure=3e5,
    temperature=300,
)

VolumeFluid = FluidLookup(
    "Volume Fluid",
    MixtureNetwork,
    SourceFluid.composition,
    pressure=2e5,
    temperature=300,
)

SeparatorOutlet1Fluid = FluidLookup(
    "Separator Outlet 1 Fluid",
    MixtureNetwork,
    SeparatorOutlet1Composition,
    pressure=VolumeFluid.pressure,
    temperature=VolumeFluid.temperature,
)

SeparatorOutlet2Fluid = FluidLookup(
    "Separator Outlet 2 Fluid",
    MixtureNetwork,
    SeparatorOutlet2Composition,
    pressure=VolumeFluid.pressure,
    temperature=VolumeFluid.temperature,
)

# ------------------------------------------------------------------
# Components after lookups
# ------------------------------------------------------------------

Inlet = DischargeCoefficient(
    "Inlet",
    MixtureNetwork,
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=VolumeFluid.pressure,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
)

Separator = FlowSplitter(
    "Separator",
    MixtureNetwork,
    pressure=VolumeFluid.pressure,
    volume=1,
    mass_flow_in=Inlet.mass_flow,
    composition=VolumeFluid.composition,
    composition_out1=SeparatorOutlet1Composition,
    composition_out2=SeparatorOutlet2Composition,
)

Outlet1 = DischargeCoefficient(
    "Outlet 1",
    MixtureNetwork,
    upstream_pressure=Separator.pressure,
    downstream_pressure=101325,
    density=SeparatorOutlet1Fluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A / 4,
    mass_flow=Separator.mass_flow_out1,
)

Outlet2 = DischargeCoefficient(
    "Outlet 2",
    MixtureNetwork,
    upstream_pressure=Separator.pressure,
    downstream_pressure=101325,
    density=SeparatorOutlet2Fluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=Separator.mass_flow_out2,
)

ArgonBalance = Balance(
    "Argon Balance",
    MixtureNetwork,
    variable=Outlet1.discharge_coefficient,
    function=Separator.composition_out2["Argon"],
)


# ------------------------------------------------------------------
# Solve
# ------------------------------------------------------------------

solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
)