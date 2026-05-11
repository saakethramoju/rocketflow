from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network



class GeneralFluidLookupfromPT(Component):
    """
    General thermodynamic property lookup component.

    Updates a persistent Fluid object from pressure and temperature
    States, then exposes all Fluid properties dynamically as derived
    State objects.

    Example
    -------
    fluid_lookup.density.value
    fluid_lookup.cp.value
    fluid_lookup.kinematic_viscosity.value
    """

    def __init__(
        self,
        name: str,
        network: Network,
        fluid: str,
        pressure: State,
        temperature: State,
    ):
        self.initialize_component(name, network)

        self.fluid = fluid
        self.pressure = pressure
        self.temperature = temperature

        # Persistent Fluid wrapper reused every evaluation.
        self._Fluid = Fluid(
            self.fluid,
            P=self.pressure.value,
            T=self.temperature.value,
        )

        # Cache of dynamically-created derived States.
        # Prevents creating duplicate State objects.
        self._property_states: dict[str, State] = {}

    def evaluate_states(self) -> None:
        """
        Update Fluid thermodynamic state from current P and T.
        """
        self._Fluid.pressure = self.pressure.value
        self._Fluid.temperature = self.temperature.value

    def __getattr__(self, name: str) -> State:
        """
        Dynamically expose Fluid properties as derived State objects.

        If the requested attribute is a property on the Fluid class,
        a derived State is lazily created and cached.

        Example
        -------
        fluid_lookup.density
        fluid_lookup.specific_heat
        fluid_lookup.dynamic_viscosity
        """

        # Avoid recursion during initialization.
        if "_Fluid" not in self.__dict__:
            raise AttributeError(name)

        # Check if attribute exists as a Fluid property.
        if hasattr(type(self._Fluid), name):
            attr = getattr(type(self._Fluid), name)

            if isinstance(attr, property):

                # Create State only once.
                if name not in self._property_states:
                    self._property_states[name] = State(
                        expr=lambda prop=name: getattr(self._Fluid, prop)
                    )

                return self._property_states[name]

        raise AttributeError(
            f"{self.__class__.__name__!s} has no attribute {name!r}"
        )

    @property
    def iteration_variables(self) -> list[State]:
        """No iteration variables."""
        return []

    @property
    def residuals(self) -> list[float]:
        """No residual equations."""
        return []



class DensityfromPT(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 temperature: State,
                 density: State | None = None):
        self.initialize_component(name, network)
        self.fluid = fluid
        self.pressure = pressure
        self.temperature = temperature

        if density is None:
            self.density = State()
        else:
            self.density = density

        self._Fluid = Fluid(self.fluid, P=self.pressure.value, T=self.temperature.value)
        self.density.value = self._Fluid.density

    def evaluate_states(self) -> None:
        self._Fluid.pressure = self.pressure.value
        self._Fluid.temperature = self.temperature.value
        self.density.value = self._Fluid.density

    @property
    def iteration_variables(self) -> list[State]:
        return []
    
    @property
    def residuals(self) -> list[float]:
        return []