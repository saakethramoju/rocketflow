from typing import Optional
from Global import SIUnits
from .Component import Component



class Branch(Component):

    def __init__(
        self,
        name: Optional[str] = None,
        mass_flow: float = 0.0,
        energy_flow: float = 0.0
    ):
        super().__init__(name)

        self.mass_flow = mass_flow
        self.energy_flow = energy_flow

        self.system_variables = {
            "mass_flow": {
                "units": SIUnits.MASS_FLOW.value,
            },
            "energy_flow": {
                "units": SIUnits.ENERGY_FLOW.value,
            }
        }



class Line(Branch):

    _counter = 0

    def __init__(
        self,
        name: Optional[str] = None,
        cross_sectional_area: Optional[float] = None,
        length: Optional[float] = None,
        mass_flow: float = 0.0,
        energy_flow: float = 0.0
    ):
        super().__init__(name, mass_flow, energy_flow)

        self.cross_sectional_area = cross_sectional_area
        self.length = length

        if name is None:
            self.name = f"Line {Line._counter}"
            Line._counter += 1

        if self.length is not None and self.length < 0:
            raise ValueError(f"Length of {self.name} must be nonnegative!")

        if self.cross_sectional_area is not None and self.cross_sectional_area <= 0:
            raise ValueError(
                f"Cross-sectional area of {self.name} must be positive!"
            )

        self.configuration_variables = {
            "length": {
                "units": SIUnits.LENGTH.value,
            },
            "cross_sectional_area": {
                "units": SIUnits.AREA.value,
            }
        }

    def __repr__(self):
        name_str = f"'{self.name}'"

        length_str = f"{self.length:.3e} m" if self.length is not None else "None"
        area_str = (
            f"{self.cross_sectional_area:.3e} m^2"
            if self.cross_sectional_area is not None else "None"
        )

        return (
            f"Line(name={name_str}, "
            f"length={length_str}, "
            f"cross_sectional_area={area_str}, "
            f"mass_flow={self.mass_flow:.3e} kg/s, "
            f"energy_flow={self.energy_flow:.3e} J/s)"
        )