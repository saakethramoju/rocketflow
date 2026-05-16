from System import *
from Solvers import *
from Utilities import Fluid

from constants import *

# --- Network Definition ---
SimpleNetwork = Network("Simple Network")

gas = IdealGasLookup(
    "Tank Gas",
    SimpleNetwork,
    "N2",
    pressure=3e5,
    temperature=300,
    flash_values=("pressure", "enthalpy"),
)

print(gas.enthalpy)