from enum import Enum

class SIUnits(Enum):
    PRESSURE = "Pa"
    TIME = "s"
    MASS = "kg"
    MASS_FLOW = "kg/s"
    LENGTH = "m"
    AREA = "m^2"
    VOLUME = "m^3"
    DENSITY = "kg/m^3"
    ENERGY = "J"
    ENERGY_FLOW = "J/s"


class UStoSI:
    PSI_TO_PA = 6894.757293168
    IN2_TO_M2 = 0.00064516
