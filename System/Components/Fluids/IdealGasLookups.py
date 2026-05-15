from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component, State
from Utilities import Fluid, IdealGas

if TYPE_CHECKING:
    from System import Network


class IdealGasLookup(Component):

    _THERMO_NAMES = (
        "pressure",
        "temperature",
        "enthalpy",
        "internal_energy",
        "density",
    )

    _PRESSURE_REQUIRED_PROPERTIES = {
        "density",
        "specific_volume",
        "entropy",
        "free_energy",
        "gibbs_energy",
    }

    _FLASH_PAIR_SETTERS = {
        frozenset(("pressure", "temperature")): ("pressure_temperature", ("pressure", "temperature")),
        frozenset(("pressure", "enthalpy")): ("pressure_enthalpy", ("pressure", "enthalpy")),
        frozenset(("pressure", "internal_energy")): ("pressure_internal_energy", ("pressure", "internal_energy")),
        frozenset(("pressure", "density")): ("pressure_density", ("pressure", "density")),
        frozenset(("density", "temperature")): ("density_temperature", ("density", "temperature")),
        frozenset(("density", "enthalpy")): ("density_enthalpy", ("density", "enthalpy")),
        frozenset(("density", "internal_energy")): ("density_internal_energy", ("density", "internal_energy")),
    }

    _SINGLE_FLASH_NAMES = {
        "temperature",
        "enthalpy",
        "internal_energy",
    }

    def __init__(
        self,
        name: str,
        network: Network,
        fluid: str,
        pressure: State | float | None = None,
        temperature: State | float | None = None,
        enthalpy: State | float | None = None,
        internal_energy: State | float | None = None,
        density: State | float | None = None,
        reference_temperature: float = 298.15,
        reference_pressure: float = 101325.0,
        reference_enthalpy: float = 0.0,
        reference_internal_energy: float = 0.0,
        **property_states: State,
    ):

        _input_map = {
            "pressure": pressure,
            "temperature": temperature,
            "enthalpy": enthalpy,
            "internal_energy": internal_energy,
            "density": density,
        }

        self.setup()

        if hasattr(self, "_input_map"):
            delattr(self, "_input_map")

        if hasattr(self, "property_states"):
            delattr(self, "property_states")

        self.reference_temperature = float(self.reference_temperature.value)
        self.reference_pressure = float(self.reference_pressure.value)
        self.reference_enthalpy = float(self.reference_enthalpy.value)
        self.reference_internal_energy = float(
            self.reference_internal_energy.value
        )

        self._reference_IdealGas = IdealGas(
            self.fluid,
            pressure=self.reference_pressure,
            temperature=self.reference_temperature,
        )

        provided_names = [
            prop_name
            for prop_name in self._THERMO_NAMES
            if _input_map[prop_name] is not None
        ]

        if len(provided_names) == 0:
            raise ValueError(
                "IdealGasLookup requires at least one thermodynamic input."
            )

        first_name = provided_names[0]

        if first_name in {"pressure", "density"}:
            if len(provided_names) < 2:
                raise ValueError(
                    f"{first_name} cannot define an ideal-gas state by itself."
                )

            self._flash_names = provided_names[:2]

        else:
            self._flash_names = [first_name]

        self._validate_flash_names()

        self._IdealGas = IdealGas(
            self.fluid,
            **self._ideal_gas_constructor_kwargs(),
        )

        self._property_states: dict[str, State] = {}
        self._external_property_names: set[str] = set()

        for output_name in provided_names[len(self._flash_names):]:
            state = getattr(self, output_name)
            self._property_states[output_name] = state
            self._external_property_names.add(output_name)

        for prop_name in self._THERMO_NAMES:
            if _input_map[prop_name] is None and hasattr(self, prop_name):
                delattr(self, prop_name)

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

            if not self._is_ideal_gas_property(prop_name):
                raise AttributeError(
                    f"{prop_name!r} is not a valid IdealGas property."
                )

            self._property_states[prop_name] = state
            self._external_property_names.add(prop_name)

    def pre_evaluation(self):
        self.evaluate_states()

    def evaluate_states(self) -> None:
        self._set_ideal_gas_from_flash()

        for prop_name in self._external_property_names:
            self._property_states[prop_name].value = self._get_property(prop_name)

    def __getattr__(self, name: str) -> State:

        if "_IdealGas" not in self.__dict__:
            raise AttributeError(name)

        if not self._is_ideal_gas_property(name):
            raise AttributeError(
                f"{self.__class__.__name__!s} has no attribute {name!r}"
            )

        if self._requires_pressure(name) and self._IdealGas.pressure is None:
            raise AttributeError(
                f"{name!r} requires pressure, but pressure is not available."
            )

        if name not in self._property_states:
            self._property_states[name] = State(
                expr=lambda prop=name: self._get_property(prop)
            )

        return self._property_states[name]

    def _validate_flash_names(self) -> None:
        if len(self._flash_names) == 1:
            if self._flash_names[0] not in self._SINGLE_FLASH_NAMES:
                raise ValueError(
                    f"Unsupported IdealGas flash input: {self._flash_names[0]!r}."
                )
            return

        key = frozenset(self._flash_names)

        if key not in self._FLASH_PAIR_SETTERS:
            raise ValueError(
                f"Unsupported IdealGas flash pair: {sorted(self._flash_names)}."
            )

    def _ideal_gas_constructor_kwargs(self) -> dict:
        return {
            flash_name: self._to_ideal_basis(
                flash_name,
                getattr(self, flash_name).value,
            )
            for flash_name in self._flash_names
        }

    def _set_ideal_gas_from_flash(self) -> None:
        if len(self._flash_names) == 1:
            flash_name = self._flash_names[0]

            setattr(
                self._IdealGas,
                flash_name,
                self._to_ideal_basis(
                    flash_name,
                    getattr(self, flash_name).value,
                ),
            )

            return

        setter_name, ordered_names = self._FLASH_PAIR_SETTERS[
            frozenset(self._flash_names)
        ]

        setattr(
            self._IdealGas,
            setter_name,
            tuple(
                self._to_ideal_basis(
                    prop_name,
                    getattr(self, prop_name).value,
                )
                for prop_name in ordered_names
            ),
        )

    def _to_ideal_basis(self, name: str, value: float) -> float:
        if name == "enthalpy":
            return (
                self._reference_IdealGas.enthalpy
                + float(value)
                - self.reference_enthalpy
            )

        if name == "internal_energy":
            return (
                self._reference_IdealGas.internal_energy
                + float(value)
                - self.reference_internal_energy
            )

        return float(value)

    def _from_ideal_basis(self, name: str, value: float) -> float:
        if name == "enthalpy":
            return (
                self.reference_enthalpy
                + float(value)
                - self._reference_IdealGas.enthalpy
            )

        if name == "internal_energy":
            return (
                self.reference_internal_energy
                + float(value)
                - self._reference_IdealGas.internal_energy
            )

        return float(value)

    def _get_property(self, name: str):
        if self._requires_pressure(name) and self._IdealGas.pressure is None:
            raise ValueError(
                f"{name!r} requires pressure, but pressure is not available."
            )

        value = getattr(self._IdealGas, name)

        return self._from_ideal_basis(name, value)

    def _requires_pressure(self, name: str) -> bool:
        return name in self._PRESSURE_REQUIRED_PROPERTIES

    def _is_ideal_gas_property(self, name: str) -> bool:
        return isinstance(
            getattr(IdealGas, name, None),
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
            "IdealGas",
            "_IdealGas",
            "reference_IdealGas",
            "_reference_IdealGas",
            "input_map",
            "_input_map",
        }
    
    

class Helium(IdealGasLookup):

    def __init__(
        self,
        name: str,
        network: Network,
        pressure: State | float | None = None,
        temperature: State | float | None = None,
        enthalpy: State | float | None = None,
        internal_energy: State | float | None = None,
        density: State | float | None = None,
        reference_temperature: float = 298.15,
        reference_pressure: float = 101325.0,
        **property_states: State,
    ):

        reference_fluid = Fluid(
            "Helium",
            pressure=reference_pressure,
            temperature=reference_temperature,
        )

        super().__init__(
            name=name,
            network=network,
            fluid="Helium",
            pressure=pressure,
            temperature=temperature,
            enthalpy=enthalpy,
            internal_energy=internal_energy,
            density=density,
            reference_temperature=reference_temperature,
            reference_pressure=reference_pressure,
            reference_enthalpy=reference_fluid.enthalpy,
            reference_internal_energy=reference_fluid.internal_energy,
            **property_states,
        )