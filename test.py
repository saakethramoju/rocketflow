from System import *
from Solvers import *

from constants import *


TubeNetwork = Network("Tube Flow")

SourceGas = IdealGasLookup(
    "Source Gas",
    TubeNetwork,
    {"N2": 0.78, "O2": 0.22},
    pressure=50 * PSIA_TO_PA,
    temperature=299.817,
)

print(SourceGas.enthalpy)
print(SourceGas.composition["Nitrogen"])
SourceGas.composition["Nitrogen"].value = 0.76
SourceGas.composition["Oxygen"].value = 1 - SourceGas.composition["Nitrogen"].value
print(SourceGas.composition)
SourceGas.evaluate_states()
print(SourceGas.enthalpy)