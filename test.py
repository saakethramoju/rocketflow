from System import *
from Solvers import *

from constants import *


TubeNetwork = Network("Tube Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    TubeNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=299.817,
)


ShockTube = StationaryNormalShock(
    "Shock Tube",
    TubeNetwork,
    specific_heat_ratio=SourceGas.specific_heat_ratio,
    #upstream_mach_number=1.4
    static_pressure_ratio=0.1
)

solution = SteadyState(TubeNetwork).solve(
    return_type="dataframe",
    verbose=True,
    print_solution=True
)