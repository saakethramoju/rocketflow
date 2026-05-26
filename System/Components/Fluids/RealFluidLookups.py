from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid, FluidRegistry

if TYPE_CHECKING:
    from System import Network

class FluidLookup(Component):

    _THERMO_NAMES = (
        "pressure",
        "temperature",
        "enthalpy",
        "quality",
        "density",
        "internal_energy",
    )

    _FLASH_PAIR_SETTERS = {
        frozenset(("pressure", "temperature")): ("pressure_temperature", ("pressure", "temperature")),
        frozenset(("pressure", "enthalpy")): ("pressure_enthalpy", ("pressure", "enthalpy")),
        frozenset(("pressure", "quality")): ("pressure_quality", ("pressure", "quality")),
        frozenset(("temperature", "quality")): ("temperature_quality", ("temperature", "quality")),
        frozenset(("density", "internal_energy")): ("density_internal_energy", ("density", "internal_energy")),
        frozenset(("pressure", "density")): ("pressure_density", ("pressure", "density")),
        frozenset(("pressure", "internal_energy")): ("pressure_internal_energy", ("pressure", "internal_energy")),
        frozenset(("temperature", "density")): ("temperature_density", ("temperature", "density")),
        frozenset(("density", "enthalpy")): ("density_enthalpy", ("density", "enthalpy")),
        frozenset(("temperature", "enthalpy")): ("temperature_enthalpy", ("temperature", "enthalpy")),
    }

    def __init__(
        self,
        name: str,
        network: Network,
        fluid: str,
        pressure: State | float | None = None,
        temperature: State | float | None = None,
        enthalpy: State | float | None = None,
        quality: State | float | None = None,
        density: State | float | None = None,
        internal_energy: State | float | None = None,
        flash_values: tuple[str, str] | None = None,
        **property_states: State,
    ):

        _input_map = {
            "pressure": pressure,
            "temperature": temperature,
            "enthalpy": enthalpy,
            "quality": quality,
            "density": density,
            "internal_energy": internal_energy,
        }

        self.setup()

        if hasattr(self, "_input_map"):
            delattr(self, "_input_map")

        if hasattr(self, "property_states"):
            delattr(self, "property_states")

        self._coolprop_fluid = FluidRegistry.coolprop_name(self.fluid)

        provided_names = [
            prop_name
            for prop_name in self._THERMO_NAMES
            if _input_map[prop_name] is not None
        ]

        if len(provided_names) < 2:
            raise ValueError(
                "FluidLookup requires at least two provided thermodynamic inputs "
                "so the initial fluid state can be defined."
            )

        # ---------- actual flash values used during evaluation ----------
        if flash_values is None:
            self._flash_names = provided_names[:2]
        else:
            if not isinstance(flash_values, tuple) or len(flash_values) != 2:
                raise ValueError(
                    "flash_values must be None or a tuple of two property names, "
                    "for example ('pressure', 'enthalpy')."
                )

            self._flash_names = list(flash_values)

            invalid_flash_values = [
                name for name in self._flash_names
                if name not in self._THERMO_NAMES
            ]

            if invalid_flash_values:
                raise ValueError(
                    f"Invalid flash_values: {invalid_flash_values}. "
                    f"Valid names are: {list(self._THERMO_NAMES)}."
                )

        flash_key = frozenset(self._flash_names)

        if flash_key not in self._FLASH_PAIR_SETTERS:
            raise ValueError(
                f"Unsupported Fluid flash pair: {sorted(self._flash_names)}. "
                f"Supported pairs are: {Fluid.available_flash_pairs()}."
            )

        # ---------- initial state used to initialize missing flash values ----------
        initial_flash_names = provided_names[:2]
        initial_flash_key = frozenset(initial_flash_names)

        if initial_flash_key not in self._FLASH_PAIR_SETTERS:
            raise ValueError(
                f"Unsupported initial Fluid flash pair: {sorted(initial_flash_names)}. "
                f"Supported pairs are: {Fluid.available_flash_pairs()}."
            )

        initial_setter_name, initial_ordered_names = self._FLASH_PAIR_SETTERS[
            initial_flash_key
        ]

        self._Fluid = Fluid(
            self._coolprop_fluid,
            **{
                name: getattr(self, name).value
                for name in initial_ordered_names
            },
        )

        self._property_states: dict[str, State] = {}
        self._external_property_names: set[str] = set()

        # ---------- flash properties are owned assignable States ----------
        for flash_name in self._flash_names:
            initial_value = getattr(self._Fluid, flash_name)

            state = getattr(self, flash_name, None)

            if hasattr(state, "is_assigned"):
                if not state.is_assigned:
                    state.value = initial_value
            else:
                setattr(self, flash_name, State(initial_value))

        # ---------- delete unprovided non-flash placeholders ----------
        # This lets __getattr__ dynamically create derived property States.
        for prop_name in self._THERMO_NAMES:
            if prop_name in self._flash_names:
                continue

            if _input_map[prop_name] is None and prop_name in self.__dict__:
                delattr(self, prop_name)

        # ---------- provided non-flash thermo states become output States ----------
        for prop_name in self._THERMO_NAMES:
            if prop_name in self._flash_names:
                continue

            if prop_name in self.__dict__:
                self._property_states[prop_name] = getattr(self, prop_name)
                self._external_property_names.add(prop_name)

        for prop_name, state in property_states.items():

            state = self.initialize_attribute(state)

            if not isinstance(state, State):
                raise TypeError(
                    f"{prop_name!r} must be a State, "
                    f"got {type(state).__name__}."
                )

            if prop_name in self._flash_names:
                raise ValueError(
                    f"{prop_name!r} is already being used as a flash input and "
                    f"cannot also be used as an output property State."
                )

            if not self._is_fluid_property(prop_name):
                raise AttributeError(
                    f"{prop_name!r} is not a valid Fluid property."
                )

            self._property_states[prop_name] = state
            self._external_property_names.add(prop_name)

    def pre_evaluation(self):
        self.evaluate_states()

    def evaluate_states(self) -> None:
        self._set_fluid_from_flash()

        for prop_name in self._external_property_names:
            self._property_states[prop_name].value = getattr(
                self._Fluid,
                prop_name,
            )

    def __getattr__(self, name: str) -> State:

        if "_Fluid" not in self.__dict__:
            raise AttributeError(name)

        if not self._is_fluid_property(name):
            raise AttributeError(
                f"{self.__class__.__name__!s} has no attribute {name!r}"
            )

        if name not in self._property_states:
            self._property_states[name] = State._derived(
                lambda prop=name: getattr(self._Fluid, prop)
            )

        return self._property_states[name]

    def _set_fluid_from_flash(self) -> None:

        setter_name, ordered_names = self._FLASH_PAIR_SETTERS[
            frozenset(self._flash_names)
        ]

        setattr(
            self._Fluid,
            setter_name,
            tuple(
                getattr(self, prop_name).value
                for prop_name in ordered_names
            ),
        )

    def _is_fluid_property(self, name: str) -> bool:
        return isinstance(
            getattr(Fluid, name, None),
            property,
        )

    @property
    def ignored_export_attributes(self) -> set[str]:
        return super().ignored_export_attributes | {
            "property_states",
            "_property_states",
            "external_property_names",
            "_external_property_names",
            "flash_names",
            "_flash_names",
            "Fluid",
            "_Fluid",
            "input_map",
            "_input_map",
            "coolprop_fluid",
            "_coolprop_fluid",
        }

class DensityfromPT(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 temperature: State,
                 density: State | None = None):
        
        self.setup()
        self._Fluid = Fluid(self.fluid, pressure=self.pressure.value, temperature=self.temperature.value)
        self.density.value = self._Fluid.density

    def evaluate_states(self) -> None:
        self._Fluid.pressure = self.pressure.value
        self._Fluid.temperature = self.temperature.value
        self.density.value = self._Fluid.density




class EnthalpyfromPT(Component):

    def __init__(self, 
                 name: str,
                 network: Network,
                 fluid: str,
                 pressure: State,
                 temperature: State,
                 enthalpy: State | None = None):
        
        self.setup()
        self._Fluid = Fluid(self.fluid, pressure=self.pressure.value, temperature=self.temperature.value)
        self.enthalpy.value = self._Fluid.density

    def evaluate_states(self) -> None:
        self._Fluid.pressure = self.pressure.value
        self._Fluid.temperature = self.temperature.value
        self.enthalpy.value = self._Fluid.density
