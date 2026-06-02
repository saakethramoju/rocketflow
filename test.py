from System import *
from Solvers import *

from constants import *
'''
MixtureNetwork = Network("Mixture Flow")

# ------------------------------------------------------------------
# Geometry
# ------------------------------------------------------------------

L = 60 * IN_TO_M
D = 3 * IN_TO_M
A = (np.pi / 4) * D**2

# ------------------------------------------------------------------
# Fluid lookups
# ------------------------------------------------------------------


SourceFluid1 = FluidLookup(
    "Source Fluid 1",
    MixtureNetwork,
    {"gn2": 1.0},
    pressure=5e5,
    temperature=300,
)


SourceFluid2 = IdealGasLookup(
    "Source Fluid 2",
    MixtureNetwork,
    {"o2": 1.0},
    pressure=5e5,
    temperature=300,
)


MixerFluid = FluidLookup(
    "Source Fluid 2",
    MixtureNetwork,
    Composition("o2"), # fix so this can be empty
    pressure=4.7e5,
    temperature=300,
    flash_values=("pressure", "enthalpy")
)

DrainFluid = FluidLookup(
    "Drain",
    MixtureNetwork,
    MixerFluid.composition,
    pressure=101325,
    temperature=290,
    flash_values=("pressure", "enthalpy")
)




# ------------------------------------------------------------------
# Components
# ------------------------------------------------------------------

Inlet1 = DischargeCoefficient(
    "Inlet 1",
    MixtureNetwork,
    upstream_pressure=SourceFluid1.pressure,
    downstream_pressure=MixerFluid.pressure,
    density=SourceFluid1.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
)


Inlet2 = DischargeCoefficient(
    "Inlet 2",
    MixtureNetwork,
    upstream_pressure=SourceFluid2.pressure,
    downstream_pressure=MixerFluid.pressure,
    density=SourceFluid2.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
)

Mixer = FlowMixer(
    "Mixer",
    MixtureNetwork,
    pressure=MixerFluid.pressure,
    volume=1,
    mass_flow_in1=Inlet1.mass_flow,
    mass_flow_in2=Inlet2.mass_flow,
    mass_flow_out=3.7,
    composition_in1=SourceFluid1.composition,
    composition_in2=SourceFluid2.composition,
    composition=MixerFluid.composition,
    total_enthalpy_in1=SourceFluid1.enthalpy,
    total_enthalpy_in2=SourceFluid2.enthalpy,
    enthalpy=MixerFluid.enthalpy,
)

"""
Outlet = DischargeCoefficient(
    "Outlet 1",
    MixtureNetwork,
    upstream_pressure=Mixer.pressure,
    downstream_pressure=101325,
    density=MixerFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=Mixer.mass_flow_out,
)
"""


Outlet = CompressibleFlowTube(
    "Outlet",
    MixtureNetwork,
    mass_flow=Mixer.mass_flow_out,
    upstream_static_pressure=Mixer.pressure,
    upstream_density=MixerFluid.density,
    upstream_static_temperature=MixerFluid.temperature,
    downstream_static_pressure=DrainFluid.pressure,
    downstream_static_temperature=DrainFluid.temperature,
    downstream_density=DrainFluid.density,
    length=L,
    inner_diameter=D,
    friction_factor=2e-5,
    upstream_static_enthalpy=MixerFluid.enthalpy,
    total_enthalpy=Mixer.total_enthalpy_out
)

# ------------------------------------------------------------------
# Solve
# ------------------------------------------------------------------

solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
    jacobian_method='2-point',
)
'''



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
    flash_values=("pressure", "enthalpy")
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
    total_enthalpy_in=SourceFluid.enthalpy,
    enthalpy=VolumeFluid.enthalpy,
    total_enthalpy_out1=SeparatorOutlet1Fluid.enthalpy,
    total_enthalpy_out2=SeparatorOutlet2Fluid.enthalpy
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
