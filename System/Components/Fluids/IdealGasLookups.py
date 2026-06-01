from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from System import Component, State, Composition
from Utilities import Fluid, IdealGas, FluidRegistry

if TYPE_CHECKING:
    from System import Network


FluidInput = str | dict[str, State | float] | Composition


class IdealGasLookup(Component):
    """
    PYroMat-backed ideal-gas property lookup component.
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
        fluid: FluidInput,
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

        initial_fluid = self.fluid
        self.composition = self._initialize_composition(initial_fluid)
        self.fluid = self.composition

        self._last_composition_values: tuple[float, ...] | None = None

        self._coolprop_fluid = None
        self._pyromat_fluid = None

        self._reference_enthalpy = None
        self._reference_internal_energy = None
        self._reference_IdealGas = None
        self._IdealGas = None

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
        self._initial_flash_names = initial_flash_names

        self._last_flash_values: tuple[float, ...] | None = None
        self._property_cache: dict[str, float] = {}

        self._property_states: dict[str, State] = {}
        self._external_property_names: set[str] = set()

        # Initialize now if possible; otherwise defer until composition is valid.
        if self.composition.is_assigned and self._composition_is_valid():
            self._initialize_backend()

        # Flash properties are owned assignable States.
        for flash_name in self._flash_names:
            state = getattr(self, flash_name, None)

            if hasattr(state, "is_assigned"):
                if self._IdealGas is not None and not state.is_assigned:
                    state.value = self._get_property(flash_name)
            else:
                setattr(self, flash_name, State(state))

        # Delete unprovided non-flash placeholders.
        for prop_name in self._THERMO_NAMES:
            if prop_name in self._flash_names:
                continue

            if _input_map[prop_name] is None and prop_name in self.__dict__:
                delattr(self, prop_name)

        # Provided non-flash thermo states become output States.
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
        if not self._ensure_backend_initialized():
            return

        self._set_ideal_gas_from_composition()
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

        if name not in self._property_states:
            self._property_states[name] = State._derived(
                lambda prop=name: self._get_property(prop)
            )

        return self._property_states[name]

    def _initialize_backend(self) -> None:
        self._validate_composition_support()

        self._coolprop_fluid = self._coolprop_argument_from_composition()
        self._pyromat_fluid = self._pyromat_argument_from_composition()

        reference_fluid = Fluid(
            self._coolprop_fluid,
            basis="mass",
            pressure=self._REFERENCE_PRESSURE,
            temperature=self._REFERENCE_TEMPERATURE,
        )

        self._reference_enthalpy = reference_fluid.enthalpy
        self._reference_internal_energy = reference_fluid.internal_energy

        self._reference_IdealGas = IdealGas(
            self._pyromat_fluid,
            basis="mass",
            pressure=self._REFERENCE_PRESSURE,
            temperature=self._REFERENCE_TEMPERATURE,
        )

        self._IdealGas = IdealGas(
            self._pyromat_fluid,
            basis="mass",
            **{
                flash_name: self._to_ideal_basis(
                    flash_name,
                    getattr(self, flash_name).value,
                )
                for flash_name in self._initial_flash_names
            },
        )

        self._last_flash_values = None
        self._property_cache.clear()

    def _ensure_backend_initialized(self) -> bool:
        if self._IdealGas is not None:
            return True

        if not self.composition.is_assigned:
            return False

        if not self._composition_is_valid():
            return False

        self._initialize_backend()
        return True

    def _composition_is_valid(self) -> bool:
        if not self.composition.is_assigned:
            return False

        values = tuple(
            self.composition[species].value
            for species in self.composition.species
        )

        return np.isclose(sum(values), 1.0, rtol=0.0, atol=1e-6)

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
        if not self._ensure_backend_initialized():
            raise ValueError(
                f"{self.name}: cannot evaluate {name!r} because the "
                "ideal-gas composition is not initialized yet."
            )

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

    def _initialize_composition(self, fluid: FluidInput) -> Composition:

        if isinstance(fluid, Composition):
            return fluid

        try:
            composition = Composition(fluid)
        except Exception as e:
            raise ValueError(
                f"{self.name}: invalid fluid input {fluid!r}. "
                "Expected a fluid name, a species-fraction dictionary, "
                "or a Composition object."
            ) from e

        if not composition.is_assigned:
            raise ValueError(
                f"{self.name}: composition must contain at least one species."
            )

        return composition

    def _validate_composition_support(self) -> None:

        for species in self.composition.species:
            if not FluidRegistry.supports_both(species):
                raise ValueError(
                    f"{species!r} must be supported by both CoolProp and PYroMat "
                    f"because IdealGasLookup uses CoolProp reference enthalpy/internal "
                    f"energy and PYroMat ideal-gas properties."
                )

    def _coolprop_argument_from_composition(self) -> str | dict[str, float]:

        values = self.composition.values

        if len(values) == 1:
            species = next(iter(values))
            return FluidRegistry.coolprop_name(species)

        return FluidRegistry.coolprop_mixture_dict(values)

    def _pyromat_argument_from_composition(self) -> str | dict[str, float]:

        values = self.composition.values

        if len(values) == 1:
            species = next(iter(values))
            return FluidRegistry.pyromat_name(species, include_prefix=True)

        return FluidRegistry.pyromat_mixture_dict(
            values,
            include_prefix=True,
        )

    def _composition_values(self) -> tuple[float, ...]:
        return tuple(
            self.composition[species].value
            for species in self.composition.species
        )

    def _composition_values_unchanged(
        self,
        composition_values: tuple[float, ...],
        rtol: float = 1e-10,
        atol: float = 1e-12,
    ) -> bool:

        if self._last_composition_values is None:
            return False

        return all(
            np.isclose(current, previous, rtol=rtol, atol=atol)
            for current, previous in zip(
                composition_values,
                self._last_composition_values,
            )
        )


    def _set_ideal_gas_from_composition(self) -> None:

        composition_values = self._composition_values()

        if self._composition_values_unchanged(composition_values):
            return

        total = sum(composition_values)

        if not np.isclose(total, 1.0, rtol=0.0, atol=1e-6):
            raise ValueError(
                f"{self.name}: composition mass fractions must sum to 1.0. "
                f"Got {total}."
            )

        new_coolprop_fluid = self._coolprop_argument_from_composition()
        new_pyromat_fluid = self._pyromat_argument_from_composition()

        # Species set changed, so rebuild the ideal-gas and reference backends.
        if (
            new_coolprop_fluid != self._coolprop_fluid
            or new_pyromat_fluid != self._pyromat_fluid
        ):
            self._initialize_backend()

            self._last_composition_values = composition_values
            self._last_flash_values = None
            self._property_cache.clear()
            return

        # Same species set, only fractions changed.
        if len(composition_values) > 1:
            self._IdealGas.mass_fractions = list(composition_values)
            self._reference_IdealGas.mass_fractions = list(composition_values)

            reference_fluid = Fluid(
                self._coolprop_fluid,
                basis="mass",
                pressure=self._REFERENCE_PRESSURE,
                temperature=self._REFERENCE_TEMPERATURE,
            )

            self._reference_enthalpy = reference_fluid.enthalpy
            self._reference_internal_energy = reference_fluid.internal_energy

        self._last_composition_values = composition_values
        self._last_flash_values = None
        self._property_cache.clear()

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
            "composition",
            "last_composition_values",
            "_last_composition_values",
            "initial_flash_names",
            "_initial_flash_names",
        }