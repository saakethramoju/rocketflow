from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

if TYPE_CHECKING:
    from System import Network


class GeneralFluidLookupfromPT(Component):
    """
    General thermodynamic property lookup component.

    A persistent Fluid object is updated from pressure and temperature
    States, then all Fluid properties are exposed dynamically as State
    objects.

    Users may optionally provide external States for selected properties:

        density=rho_state
        dynamic_viscosity=mu_state

    Otherwise, properties are lazily created as derived States when first
    accessed.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        fluid: str,
        pressure: State,
        temperature: State,
        **property_states: State,
    ):
        self.setup()

        # Persistent Fluid wrapper reused every evaluation.
        # Avoids repeatedly constructing Fluid objects.
        self._Fluid = Fluid(
            self.fluid,
            P=self.pressure.value,
            T=self.temperature.value,
        )

        # Stores ALL exposed property States:
        # - user-supplied States
        # - lazily-created derived States
        self._property_states: dict[str, State] = {}

        # Tracks only user-supplied States so we know which
        # ones should receive direct value assignments.
        self._external_property_names: set[str] = set()

        # Register any user-provided property States.
        for prop_name, state in property_states.items():

            if not isinstance(state, State):
                raise TypeError(
                    f"{prop_name!r} must be a State, "
                    f"got {type(state).__name__}."
                )

            # Ensure property exists on Fluid class.
            if not self._is_fluid_property(prop_name):
                raise AttributeError(
                    f"{prop_name!r} is not a valid Fluid property."
                )

            self._property_states[prop_name] = state
            self._external_property_names.add(prop_name)

    
    def pre_evaluation(self):
        """
        Included for the small chance that an external, user-defined State
        is passed in for a fluid property in this class, like density, and 
        that State's .value is used in a component before the fluid property
        could be evaluated, like if the fluid lookup was AFTER the component
        that calls .value on the external State.
        """
        self.evaluate_states()
    

    def evaluate_states(self) -> None:
        """
        Update Fluid thermodynamic state and synchronize any
        user-supplied property States.
        """

        # Update internal Fluid wrapper from current P/T.
        self._Fluid.pressure = self.pressure.value
        self._Fluid.temperature = self.temperature.value

        # Explicitly update external States supplied by user.
        #
        # Derived States do not need manual updates because
        # their expr callables evaluate dynamically.
        for prop_name in self._external_property_names:
            self._property_states[prop_name].value = getattr(
                self._Fluid,
                prop_name,
            )

    def __getattr__(self, name: str) -> State:
        """
        Dynamically expose Fluid properties as cached State objects.

        Example
        -------
        fluid_lookup.density
        fluid_lookup.specific_heat
        fluid_lookup.kinematic_viscosity
        """

        # Prevent recursion during object initialization.
        if "_Fluid" not in self.__dict__:
            raise AttributeError(name)

        # Reject invalid Fluid properties.
        if not self._is_fluid_property(name):
            raise AttributeError(
                f"{self.__class__.__name__!s} "
                f"has no attribute {name!r}"
            )

        # Lazily create derived State only once.
        #
        # Future accesses return the same exact State object.
        if name not in self._property_states:

            self._property_states[name] = State(
                expr=lambda prop=name: getattr(self._Fluid, prop)
            )

        return self._property_states[name]

    def _is_fluid_property(self, name: str) -> bool:
        """
        Return True if `name` is a property defined on Fluid.
        """

        if "_Fluid" not in self.__dict__:
            return False

        return isinstance(
            getattr(type(self._Fluid), name, None),
            property,
        )



class DensityfromPT(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 temperature: State,
                 density: State | None = None):
        
        self.setup()
        self._Fluid = Fluid(self.fluid, P=self.pressure.value, T=self.temperature.value)
        self.density.value = self._Fluid.density

    def evaluate_states(self) -> None:
        self._Fluid.pressure = self.pressure.value
        self._Fluid.temperature = self.temperature.value
        self.density.value = self._Fluid.density
