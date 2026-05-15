from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid

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

        provided_names = [
            prop_name
            for prop_name in self._THERMO_NAMES
            if _input_map[prop_name] is not None
        ]

        if len(provided_names) < 2:
            raise ValueError(
                "FluidLookup requires at least two provided thermodynamic inputs."
            )

        self._flash_names = provided_names[:2]

        flash_key = frozenset(self._flash_names)

        if flash_key not in self._FLASH_PAIR_SETTERS:
            raise ValueError(
                f"Unsupported Fluid flash pair: {sorted(self._flash_names)}. "
                f"Supported pairs are: {Fluid.available_flash_pairs()}."
            )

        self._property_states: dict[str, State] = {}
        self._external_property_names: set[str] = set()

        for output_name in provided_names[2:]:
            state = getattr(self, output_name)
            self._property_states[output_name] = state
            self._external_property_names.add(output_name)

        for prop_name in self._THERMO_NAMES:
            if _input_map[prop_name] is None and hasattr(self, prop_name):
                delattr(self, prop_name)

        self._Fluid = Fluid(
            self.fluid,
            **{
                flash_name: getattr(self, flash_name).value
                for flash_name in self._flash_names
            },
        )

        for prop_name, state in property_states.items():

            state = self.initialize_attribute(state)

            if not isinstance(state, State):
                raise TypeError(
                    f"{prop_name!r} must be a State, "
                    f"got {type(state).__name__}."
                )

            if prop_name in self._flash_names:
                raise ValueError(
                    f"{prop_name!r} was provided as a flash input and "
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
            self._property_states[name] = State(
                expr=lambda prop=name: getattr(self._Fluid, prop)
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
