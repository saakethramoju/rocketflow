from typing import Optional
from enum import Enum


class BoundaryType(Enum):
    PRESSURE = "pressure"
    MASS_FLOW = "mass_flow"


class Boundary:
    def __init__(self,
                 name: Optional[str] = None,
                 boundary_type: Optional[BoundaryType] = None,
                 pressure: Optional[float] = None,
                 mass_flow: Optional[float] = None):

        self.name = name
        self.pressure = pressure
        self.mass_flow = mass_flow

        if boundary_type is None:
            if pressure is not None:
                self.boundary_type = BoundaryType.PRESSURE
            elif mass_flow is not None:
                self.boundary_type = BoundaryType.MASS_FLOW
            else:
                self.boundary_type = BoundaryType.PRESSURE
        else:
            self.boundary_type = boundary_type

        if self.boundary_type == BoundaryType.PRESSURE and self.pressure is None:
            raise ValueError("Pressure boundary requires a pressure to be defined")

        if self.boundary_type == BoundaryType.MASS_FLOW and self.mass_flow is None:
            raise ValueError("Mass flow boundary requires a mass flow rate to be defined")

    def __repr__(self):
        cls_name = self.__class__.__name__

        if self.boundary_type == BoundaryType.PRESSURE:
            return (f"{cls_name}(name='{self.name}', "
                    f"type='pressure', "
                    f"pressure={self.pressure:.3e} Pa, "
                    f"mass_flow={self.mass_flow}) kg/s")
        else:
            return (f"{cls_name}(name='{self.name}', "
                    f"type='mass_flow', "
                    f"mass_flow={self.mass_flow:.3e} kg/s, "
                    f"pressure={self.pressure}) Pa")


class Source(Boundary):
    pressure_counter = 0
    mass_flow_counter = 0

    def __init__(self,
                 name: Optional[str] = None,
                 boundary_type: Optional[BoundaryType] = None,
                 pressure: Optional[float] = None,
                 mass_flow: Optional[float] = None):

        super().__init__(name, boundary_type, pressure, mass_flow)

        if self.name is None:
            if self.boundary_type == BoundaryType.PRESSURE:
                self.name = f"Pressure Source {Source.pressure_counter}"
                Source.pressure_counter += 1
            else:
                self.name = f"Mass Flow Source {Source.mass_flow_counter}"
                Source.mass_flow_counter += 1


class Drain(Boundary):
    pressure_counter = 0
    mass_flow_counter = 0

    def __init__(self,
                 name: Optional[str] = None,
                 boundary_type: Optional[BoundaryType] = None,
                 pressure: Optional[float] = None,
                 mass_flow: Optional[float] = None):

        super().__init__(name, boundary_type, pressure, mass_flow)

        if self.name is None:
            if self.boundary_type == BoundaryType.PRESSURE:
                self.name = f"Pressure Drain {Drain.pressure_counter}"
                Drain.pressure_counter += 1
            else:
                self.name = f"Mass Flow Drain {Drain.mass_flow_counter}"
                Drain.mass_flow_counter += 1