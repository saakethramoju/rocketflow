from typing import Optional

class Branch:
    def __init__(self,
                 name: Optional[str] = None,
                 mass_flow: float = 0.0,
                 energy_flow: float = 0.0):

        self.name = name
        self.mass_flow = mass_flow
        self.energy_flow = energy_flow

    def __repr__(self):
        name_str = f"'{self.name}'" if self.name is not None else "None"

        return (f"Branch(name={name_str}, "
                f"mass_flow={self.mass_flow:.3e} kg/s, "
                f"energy_flow={self.energy_flow:.3e}) J/s")
    

class Line(Branch):

    counter = 0

    def __init__(self, 
                 name = None, 
                 cross_sectional_area : Optional[float] = None,
                 length : Optional[float] = None,
                 mass_flow = 0, 
                 energy_flow = 0):
        super().__init__(name, mass_flow, energy_flow)
        self.cross_sectional_area = cross_sectional_area
        self.length = length

        if self.name is None:
            self.name = f"Line {Line.counter}"
            Line.counter += 1
        
        if self.length <= 0:
            raise ValueError(f"Length of {self.name} must be positive!")
        if self.cross_sectional_area <= 0:
            raise ValueError(f"Cross-sectional area of {self.name} must be positive!")

    def __repr__(self):
        name_str = f"'{self.name}'"

        return (
            f"Line(name={name_str}, "
            f"length={self.length:.3e} m, "
            f"cross-sectional area={self.cross_sectional_area:.3e} m^2, "
            f"mass_flow={self.mass_flow:.3e} kg/s, "
            f"energy_flow={self.energy_flow:.3e} J/s)"
        )