from System import *
from Solvers import *

from constants import *


MixtureNetwork = Network("Mixture Flow")

SourceFluid = IdealGasLookup(
    "Source Fluid",
    MixtureNetwork,
    {"gn2": 0.75, "O2": 0.25},
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

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

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
    composition_out1=Composition("N2"),
)

Outlet1Fluid = FluidLookup(
    "Outlet 1 Fluid",
    MixtureNetwork,
    Separator.composition_out1,
    pressure=Separator.pressure,
    temperature=VolumeFluid.temperature,
)

Outlet2Fluid = FluidLookup(
    "Outlet 2 Fluid",
    MixtureNetwork,
    Separator.composition_out2,
    pressure=Separator.pressure,
    temperature=VolumeFluid.temperature,
)

Outlet1 = DischargeCoefficient(
    "Outlet 1",
    MixtureNetwork,
    upstream_pressure=Separator.pressure,
    downstream_pressure=101325,
    density=Outlet1Fluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=Separator.mass_flow_out1,
)

Outlet2 = DischargeCoefficient(
    "Outlet 2",
    MixtureNetwork,
    upstream_pressure=Separator.pressure,
    downstream_pressure=101325,
    density=Outlet2Fluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=Separator.mass_flow_out2,
)

solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
    filename="test.xlsx",
)