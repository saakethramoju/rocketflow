from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from System import Component, State
from Utilities import Fluid, IdealGas, FluidRegistry

if TYPE_CHECKING:
    from System import Network


class IdealGasLookup(Component):
    """
    PYroMat-backed ideal-gas property lookup component.

    IdealGasLookup owns a persistent IdealGas object and updates it from either a
    single ideal-gas property, such as temperature or enthalpy, or a two-property
    flash pair, such as pressure-temperature or pressure-enthalpy. The flash inputs
    are State objects and may be solver iteration variables or shared states.

    The component uses CoolProp only once during initialization to align enthalpy
    and internal-energy reference values with the real-fluid Fluid wrapper. Runtime
    property evaluation is handled by the PYroMat-backed IdealGas object.

    To reduce repeated PYroMat work, the lookup skips unchanged flashes and caches
    property values after each flash. Cached properties are reused until the flash
    inputs change, at which point the cache is cleared and rebuilt lazily.
    """

    _REFERENCE_TEMPERATURE = 298.15
    _REFERENCE_PRESSURE = 101325.0

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
        flash_values: tuple[str, ...] | None = None,
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

        if not FluidRegistry.supports_both(self.fluid):
            raise ValueError(
                f"{self.fluid!r} must be supported by both CoolProp and "
                f"PYroMat because IdealGasLookup uses CoolProp reference "
                f"enthalpy/internal energy and PYroMat ideal-gas properties."
            )

        self._coolprop_fluid = FluidRegistry.coolprop_name(self.fluid)
        self._pyromat_fluid = FluidRegistry.pyromat_name(self.fluid)

        reference_fluid = Fluid(
            self._coolprop_fluid,
            pressure=self._REFERENCE_PRESSURE,
            temperature=self._REFERENCE_TEMPERATURE,
        )

        self._reference_enthalpy = reference_fluid.enthalpy
        self._reference_internal_energy = reference_fluid.internal_energy

        self._reference_IdealGas = IdealGas(
            self._pyromat_fluid,
            pressure=self._REFERENCE_PRESSURE,
            temperature=self._REFERENCE_TEMPERATURE,
        )

        provided_names = [
            prop_name
            for prop_name in self._THERMO_NAMES
            if _input_map[prop_name] is not None
        ]

        if len(provided_names) == 0:
            raise ValueError(
                "IdealGasLookup requires at least one thermodynamic input "
                "so the initial ideal-gas state can be defined."
            )

        if flash_values is None:
            if "pressure" in provided_names:
                if len(provided_names) < 2:
                    raise ValueError(
                        "pressure cannot define an ideal-gas state by itself."
                    )

                self._flash_names = provided_names[:2]

                if "pressure" not in self._flash_names:
                    raise ValueError(
                        "If pressure is provided, it must be one of the first "
                        "two thermodynamic inputs."
                    )

            elif "density" in provided_names:
                if len(provided_names) < 2:
                    raise ValueError(
                        "density cannot define an ideal-gas state by itself."
                    )

                self._flash_names = provided_names[:2]

                if "density" not in self._flash_names:
                    raise ValueError(
                        "If density is provided, it must be one of the first "
                        "two thermodynamic inputs."
                    )

            else:
                self._flash_names = [provided_names[0]]

        else:
            if not isinstance(flash_values, tuple) or len(flash_values) not in {1, 2}:
                raise ValueError(
                    "flash_values must be None, a tuple with one property name, "
                    "or a tuple with two property names. Examples: "
                    "('temperature',) or ('pressure', 'enthalpy')."
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

        self._validate_flash_names()

        if len(provided_names) == 1:
            initial_flash_names = [provided_names[0]]
        else:
            initial_flash_names = provided_names[:2]

        self._validate_initial_flash_names(initial_flash_names)

        self._IdealGas = IdealGas(
            self._pyromat_fluid,
            **{
                flash_name: self._to_ideal_basis(
                    flash_name,
                    getattr(self, flash_name).value,
                )
                for flash_name in initial_flash_names
            },
        )

        self._last_flash_values: tuple[float, ...] | None = None
        self._property_cache: dict[str, float] = {}

        self._property_states: dict[str, State] = {}
        self._external_property_names: set[str] = set()

        # ---------- flash properties are owned assignable States ----------
        for flash_name in self._flash_names:
            initial_value = self._get_property(flash_name)

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
            self._property_states[name] = State._derived(
                lambda prop=name: self._get_property(prop)
            )

        return self._property_states[name]

    def _validate_initial_flash_names(self, flash_names: list[str]) -> None:
        if len(flash_names) == 1:
            if flash_names[0] not in self._SINGLE_FLASH_NAMES:
                raise ValueError(
                    f"Unsupported initial IdealGas flash input: {flash_names[0]!r}."
                )
            return

        key = frozenset(flash_names)

        if key not in self._FLASH_PAIR_SETTERS:
            raise ValueError(
                f"Unsupported initial IdealGas flash pair: {sorted(flash_names)}."
            )

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

    def _set_ideal_gas_from_flash(self) -> None:
        if len(self._flash_names) == 1:
            flash_name = self._flash_names[0]

            flash_values = (
                self._to_ideal_basis(
                    flash_name,
                    getattr(self, flash_name).value,
                ),
            )

            if self._flash_values_unchanged(flash_values):
                return

            setattr(
                self._IdealGas,
                flash_name,
                flash_values[0],
            )

            self._last_flash_values = flash_values
            self._property_cache.clear()

            return

        setter_name, ordered_names = self._FLASH_PAIR_SETTERS[
            frozenset(self._flash_names)
        ]

        flash_values = tuple(
            self._to_ideal_basis(
                prop_name,
                getattr(self, prop_name).value,
            )
            for prop_name in ordered_names
        )

        if self._flash_values_unchanged(flash_values):
            return

        setattr(
            self._IdealGas,
            setter_name,
            flash_values,
        )

        self._last_flash_values = flash_values
        self._property_cache.clear()

    def _flash_values_unchanged(
        self,
        flash_values: tuple[float, ...],
        rtol: float = 1e-10,
        atol: float = 1e-12,
    ) -> bool:

        if self._last_flash_values is None:
            return False

        return all(
            np.isclose(current, previous, rtol=rtol, atol=atol)
            for current, previous in zip(
                flash_values,
                self._last_flash_values,
            )
        )

    def _to_ideal_basis(self, name: str, value: float) -> float:
        if name == "enthalpy":
            return (
                self._reference_IdealGas.enthalpy
                + float(value)
                - self._reference_enthalpy
            )

        if name == "internal_energy":
            return (
                self._reference_IdealGas.internal_energy
                + float(value)
                - self._reference_internal_energy
            )

        return float(value)

    def _from_ideal_basis(self, name: str, value: float) -> float:
        if name == "enthalpy":
            return (
                self._reference_enthalpy
                + float(value)
                - self._reference_IdealGas.enthalpy
            )

        if name == "internal_energy":
            return (
                self._reference_internal_energy
                + float(value)
                - self._reference_IdealGas.internal_energy
            )

        return float(value)

    def _get_property(self, name: str):
        if self._requires_pressure(name) and self._IdealGas.pressure is None:
            raise ValueError(
                f"{name!r} requires pressure, but pressure is not available."
            )

        if name not in self._property_cache:
            value = getattr(self._IdealGas, name)
            self._property_cache[name] = self._from_ideal_basis(name, value)

        return self._property_cache[name]

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
            "coolprop_fluid",
            "_coolprop_fluid",
            "pyromat_fluid",
            "_pyromat_fluid",
            "reference_enthalpy",
            "_reference_enthalpy",
            "reference_internal_energy",
            "_reference_internal_energy",
            "input_map",
            "_input_map",
            "last_flash_values",
            "_last_flash_values",
            "property_cache",
            "_property_cache",
        }