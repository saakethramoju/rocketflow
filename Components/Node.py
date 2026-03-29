from typing import Optional
from Global import SIUnits
from.Component import Component


class Node(Component):

    def __init__(
        self,
        name: Optional[str] = None,
        pressure: Optional[float] = None
    ):
        super().__init__(name)

        self.pressure = pressure

        if self.pressure is None:
            raise ValueError(f"{self.name} pressure / initial guess must be defined!")
        if self.pressure < 0:
            raise ValueError(f"{self.name} pressure must be nonnegative!")

        self.iteration_variables = {
            "pressure": {
                "units": SIUnits.PRESSURE.value,
            }
        }




class Volume(Node):

    _counter = 0

    def __init__(
        self,
        name: Optional[str] = None,
        pressure: Optional[float] = None,
        volume: Optional[float] = None
    ):
        super().__init__(name, pressure)

        self.volume = volume

        if name is None:
            self.name = f"Volume {Volume._counter}"
            Volume._counter += 1

        if self.volume is not None and self.volume < 0:
            raise ValueError(f"{self.name} volume must be nonnegative!")

        self.configuration_variables = {
            "volume": {
                "units": SIUnits.VOLUME.value,
            }
        }

    def __repr__(self):
        name_str = f"'{self.name}'"
        volume_str = f"{self.volume:.3e} m^3" if self.volume is not None else "None"
        return (
            f"Volume(name={name_str}, "
            f"pressure={self.pressure:.3e} Pa, "
            f"volume={volume_str})"
        )
