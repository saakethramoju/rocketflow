from System import *
from Solvers import *

from constants import *


TubeNetwork = Network("Tube Flow")

SourceGas = FluidLookup(
    "Source Gas",
    TubeNetwork,
    "gn2",
    pressure=50 * PSIA_TO_PA,
    temperature=299.817,
)

print(SourceGas.enthalpy)