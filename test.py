from System import *
from Solvers import *

from constants import *


MixtureNetwork = Network("Mixture Flow")

# ------------------------------------------------------------------
# Fluid properties
# ------------------------------------------------------------------

SourceFluid = FluidLookup(
    "Source Fluid",
    MixtureNetwork,
    #{"gn2": 0.75, "O2": 0.24, "Ar": 0.01},
    {"gn2": 0.9, "O2": 0.1},
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
    "N2",
    #Composition("Ar"), # also doesn't work idk why, some operand error
    pressure=VolumeFluid.pressure,
    temperature=VolumeFluid.temperature,
)

SeparatorOutlet2Fluid = FluidLookup(
    "Separator Outlet 2 Fluid",
    MixtureNetwork,
    {"gn2": 1.0, "O2": 0},
    #"o2", # Doesn't Work, figure out why, also 
    pressure=VolumeFluid.pressure,
    temperature=VolumeFluid.temperature,
)

# ------------------------------------------------------------------
# Geometry
# ------------------------------------------------------------------

D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

# ------------------------------------------------------------------
# Components
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
    composition_out1=SeparatorOutlet1Fluid.composition,
    composition_out2=SeparatorOutlet2Fluid.composition,
)

Outlet1 = DischargeCoefficient(
    "Outlet 1",
    MixtureNetwork,
    upstream_pressure=Separator.pressure,
    downstream_pressure=101325,
    density=SeparatorOutlet1Fluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
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

# DOESN'T WORK, FIGURE OUT WHY
'''
SeparatorBalance = Balance(
    "Pure Separator",
    MixtureNetwork,
    variable=Outlet1.discharge_coefficient,
    function=Separator.composition_out2["o2"] - 0.5
)
'''

# DOESN'T WORK, FIGURE OUT WHY
'''
SeparatorBalance = Balance(
    "Pure Separator",
    MixtureNetwork,
    variable=SourceFluid.composition["o2"],
    function=Separator.composition_out2["o2"] - 1.0
)
'''

# This seems to be the only one that works
'''
SeparatorBalance = Balance(
    "Pure Separator",
    MixtureNetwork,
    variable=SourceFluid.composition["o2"],
    function=Separator.composition_out2["o2"] - 0.5
)
'''

# ------------------------------------------------------------------
# Solve
# ------------------------------------------------------------------

solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
    filename="test.xlsx",
)



'''from System import *
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
)'''