from System import *
from Solvers import *

from constants import *


MixtureNetwork = Network("Mixture Flow")

SourceFluid = IdealGasLookup(
    "Source Fluid",
    MixtureNetwork,
    {"gn2": 0.75, "O2": 0.25},
    #"o2",
    pressure=3e5,
    temperature=300,
    #quality=0.5
)

VolumeFluid = FluidLookup(
    "Volume Fluid",
    MixtureNetwork,
    {'N2': 0.5, "gox": 0.4, "Argon": 0.1},
    pressure=2e5,
    temperature=300
)

D = 3 * IN_TO_M
A = (np.pi/4) * D**2

Inlet = DischargeCoefficient(
    "Inlet",
    MixtureNetwork,
    upstream_pressure=SourceFluid.pressure,
    downstream_pressure=VolumeFluid.pressure,
    density=SourceFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A
)

outlet_mdot1 = State()
outlet_mdot2 = State()

Vol = SimpleVolume(
    "Volume",
    MixtureNetwork,
    pressure=VolumeFluid.pressure,
    volume=1,
    mass_flow_in=Inlet.mass_flow,
    composition_in=SourceFluid.composition,
    composition=VolumeFluid.composition,
    mass_flow_out=outlet_mdot1 + outlet_mdot2
)


Outlet1 = DischargeCoefficient(
    "Outlet",
    MixtureNetwork,
    upstream_pressure=Vol.pressure,
    downstream_pressure=101325,
    density=VolumeFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A/2,
    mass_flow=outlet_mdot1
)

Outlet2 = DischargeCoefficient(
    "Outlet",
    MixtureNetwork,
    upstream_pressure=Vol.pressure,
    downstream_pressure=101325,
    density=VolumeFluid.density,
    discharge_coefficient=1,
    cross_sectional_area=A,
    mass_flow=outlet_mdot2
)

solution = SteadyState(MixtureNetwork).solve(
    return_type="dataframe",
    verbose=True,
    static=False,
    print_solution=True,
    filename='test.xlsx'
)