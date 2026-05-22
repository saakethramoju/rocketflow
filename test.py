from System import *
from Solvers import *

from constants import *


# --- Network Definition ---
PumpNetwork = Network("Pumped Sytem")

pressurant = 'methane'

FuelUllageGas = IdealGasLookup(
    "Fuel Tank Ullage Gas",
    PumpNetwork,
    pressurant,
    pressure=80 * PSIA_TO_PA,
    temperature=300,
)

print(FuelUllageGas.dynamic_viscosity.value)