from System import *
from Solvers import *

from constants import *


MixtureNetwork = Network("Mixture Flow")

SourceFluid = FluidLookup(
    "Source Fluid",
    MixtureNetwork,
    {"gn2": 0.75, "O2": 0.25},
    #"o2",
    pressure=3e5,
    temperature=300,
)

SplitterFluid = IdealGasLookup(
    "Splitter Fluid",
    MixtureNetwork,
    {"gn2": 0.75, "O2": 0.25},
    #"o2",
    pressure=2e5,
    temperature=300,
    flash_values=("pressure", "enthalpy")
)

D = 6 * IN_TO_M
A = (np.pi/4) * D**2


Line1 = DischargeCoefficient(
    "Line 1",
    MixtureNetwork,
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=SplitterFluid.pressure,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A
)

Splitter = FlowSplitter(
    "Splitter",
    MixtureNetwork,
    pressure=SplitterFluid.pressure,
    enthalpy=SplitterFluid.enthalpy,
    volume=1,
    total_enthalpy_in=SourceFluid.enthalpy,
    mass_flow_in=Line1.mass_flow,
    composition=SplitterFluid.composition,
    composition_in=SourceFluid.composition,
)


Line21 = DischargeCoefficient(
    "Line 21",
    MixtureNetwork,
    upstream_pressure=SplitterFluid.pressure,
    downstream_pressure=101325,
    density=SplitterFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=Splitter.mass_flow_out1
)

Line22 = DischargeCoefficient(
    "Line 22",
    MixtureNetwork,
    upstream_pressure=SplitterFluid.pressure,
    downstream_pressure=101325,
    density=SplitterFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=Splitter.mass_flow_out2
)

solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
    filename='test.xlsx'
)
