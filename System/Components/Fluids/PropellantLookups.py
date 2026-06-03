from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from System import Component, State, Composition
from thermoprop import Propellant

from Exceptions import InvalidThermoStateError

if TYPE_CHECKING:
    from System import Network


PropellantInput = str | Composition


class PropellantLookup(Component):
    """
    RocketProps-backed liquid propellant property lookup component.

    This lookup intentionally follows the same general pattern as FluidLookup
    and IdealGasLookup, but it is simpler because Propellant is not a full
    thermodynamic flash solver.

    Supported state inputs
    ----------------------
    Propellant only requires temperature to evaluate saturated-liquid
    properties:

        PropellantLookup(..., temperature=...)

    Pressure is optional. If pressure is provided, Propellant uses
    compressed-liquid correlations where RocketProps supports them:

        PropellantLookup(..., pressure=..., temperature=...)

    Notes
    -----
    No reference adjustment is needed because the Propellant wrapper does not
    expose reference-dependent thermodynamic properties such as enthalpy,
    internal energy, or entropy.

    Only pressure and temperature can be iteration/input states. There is no
    flash_values argument because RocketProps does not support alternate
    thermodynamic flash pairs.
    """

    _THERMO_NAMES = (
        "pressure",
        "temperature",
    )

    def __init__(
        self,
        name: str,
        network: Network,
        propellant: PropellantInput,
        temperature: State | float | None = None,
        pressure: State | float | None = None,
        **property_states: State,
    ):

        _input_map = {
            "pressure": pressure,
            "temperature": temperature,
        }

        self.setup()

        if hasattr(self, "_input_map"):
            delattr(self, "_input_map")

        if hasattr(self, "property_states"):
            delattr(self, "property_states")

        initial_propellant = self.propellant
        self.composition = self._initialize_composition(initial_propellant)
        self.propellant = self.composition

        self._rocketprops_propellant = None
        self._Propellant = None
        self._last_composition_values: tuple[float, ...] | None = None

        provided_names = [
            prop_name
            for prop_name in self._THERMO_NAMES
            if _input_map[prop_name] is not None
        ]

        if "temperature" not in provided_names:
            raise ValueError(
                "PropellantLookup requires temperature so the initial "
                "liquid propellant state can be defined."
            )

        # Propellant can be evaluated from temperature alone.
        # Pressure is optional and only enables compressed-liquid correlations.
        if "pressure" in provided_names:
            self._flash_names = ["pressure", "temperature"]
        else:
            self._flash_names = ["temperature"]

        self._last_flash_values: tuple[float, ...] | None = None
        self._property_cache: dict[str, float] = {}

        self._property_states: dict[str, State] = {}
        self._external_property_names: set[str] = set()

        # Initialize now if possible; otherwise defer until composition is valid.
        if self.composition.is_assigned and self._composition_is_valid():
            self._initialize_backend()

        # Flash states should be assignable States.
        for flash_name in self._flash_names:
            state = getattr(self, flash_name, None)

            if hasattr(state, "is_assigned"):
                if self._Propellant is not None and not state.is_assigned:
                    state.value = self._get_cached_property(flash_name)
            else:
                setattr(self, flash_name, State(state))

        # Remove unused placeholder thermo states.
        for prop_name in self._THERMO_NAMES:
            if prop_name in self._flash_names:
                continue

            if _input_map[prop_name] is None and prop_name in self.__dict__:
                delattr(self, prop_name)

        # Provided non-flash thermo states become output states.
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
                    f"{prop_name!r} is already being used as a propellant "
                    f"input and cannot also be used as an output property State."
                )

            if not self._is_propellant_property(prop_name):
                raise AttributeError(
                    f"{prop_name!r} is not a valid Propellant property."
                )

            self._property_states[prop_name] = state
            self._external_property_names.add(prop_name)

    def pre_evaluation(self):
        self.evaluate_states()

    def evaluate_states(self) -> None:
        if not self._ensure_backend_initialized():
            return

        self._set_propellant_from_composition()

        try:
            self._set_propellant_from_flash()

        except Exception as e:
            flash_state = {
                name: getattr(self, name).value
                for name in self._flash_names
            }

            composition = {
                species: state.value
                for species, state in self.composition
            }

            raise InvalidThermoStateError(
                f"{self.name}: invalid propellant state.\n"
                f"Flash variables: {flash_state}\n"
                f"Composition: {composition}"
            ) from e

        for prop_name in self._external_property_names:
            self._property_states[prop_name].value = self._get_cached_property(
                prop_name
            )

    def __getattr__(self, name: str) -> State:

        if "_Propellant" not in self.__dict__:
            raise AttributeError(name)

        if not self._is_propellant_property(name):
            raise AttributeError(
                f"{self.__class__.__name__!s} has no attribute {name!r}"
            )

        if name not in self._property_states:
            self._property_states[name] = State._derived(
                lambda prop=name: self._get_cached_property(prop)
            )

        return self._property_states[name]

    def _initialize_backend(self) -> None:
        """Create the RocketProps-backed Propellant object."""
        self._rocketprops_propellant = self._propellant_argument_from_composition()

        self._Propellant = Propellant(
            self._rocketprops_propellant,
            **{
                name: getattr(self, name).value
                for name in self._flash_names
            },
        )

        self._last_flash_values = None
        self._property_cache.clear()

    def _ensure_backend_initialized(self) -> bool:
        if self._Propellant is not None:
            return True

        if not self.composition.is_assigned:
            return False

        if not self._composition_is_valid():
            return False

        self._initialize_backend()
        return True

    def _composition_is_valid(self) -> bool:
        """
        Return True if the propellant input is valid.

        PropellantLookup only supports one RocketProps propellant name at a
        time. Named RocketProps mixtures such as MON25 or A50 should be passed
        as a single propellant name, not as a multi-species Composition.
        """
        if not self.composition.is_assigned:
            return False

        values = tuple(
            self.composition[species].value
            for species in self.composition.species
        )

        return len(values) == 1 and np.isclose(
            values[0],
            1.0,
            rtol=0.0,
            atol=1e-6,
        )

    def _set_propellant_from_flash(self) -> None:

        flash_values = tuple(
            getattr(self, prop_name).value
            for prop_name in self._flash_names
        )

        if self._flash_values_unchanged(flash_values):
            return

        if self._flash_names == ["temperature"]:
            self._Propellant.temperature = flash_values[0]

        else:
            self._Propellant.pressure_temperature = flash_values

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

    def _get_cached_property(self, name: str):

        if not self._ensure_backend_initialized():
            raise ValueError(
                f"{self.name}: cannot evaluate {name!r} because the "
                "propellant composition is not initialized yet."
            )

        if name not in self._property_cache:
            self._property_cache[name] = getattr(self._Propellant, name)

        return self._property_cache[name]

    def _is_propellant_property(self, name: str) -> bool:
        return isinstance(
            getattr(Propellant, name, None),
            property,
        )

    def _initialize_composition(self, propellant: PropellantInput) -> Composition:

        if isinstance(propellant, Composition):
            composition = propellant
        else:
            try:
                composition = Composition(propellant)
            except Exception as e:
                raise ValueError(
                    f"{self.name}: invalid propellant input {propellant!r}. "
                    "Expected a RocketProps propellant name or a single-species "
                    "Composition object."
                ) from e

        if not composition.is_assigned:
            raise ValueError(
                f"{self.name}: propellant composition must contain one species."
            )

        if len(composition.species) != 1:
            raise ValueError(
                f"{self.name}: PropellantLookup only supports one RocketProps "
                "propellant name at a time. Named mixtures such as MON25, A50, "
                "or MHF3 should be passed as a single propellant name, not as a "
                "multi-species Composition."
            )

        return composition

    def _propellant_argument_from_composition(self) -> str:

        values = self.composition.values

        if len(values) != 1:
            raise ValueError(
                f"{self.name}: PropellantLookup requires exactly one "
                "RocketProps propellant name."
            )

        propellant, fraction = next(iter(values.items()))

        if not np.isclose(fraction, 1.0, rtol=0.0, atol=1e-6):
            raise ValueError(
                f"{self.name}: single propellant fraction must be 1.0. "
                f"Got {fraction}."
            )

        return propellant

    def _composition_values(self) -> tuple[float, ...]:
        return tuple(
            self.composition[species].value
            for species in self.composition.species
        )

    def _set_propellant_from_composition(self) -> None:

        composition_values = self._composition_values()

        if self._composition_values_unchanged(composition_values):
            return

        if len(composition_values) != 1:
            raise ValueError(
                f"{self.name}: PropellantLookup only supports one RocketProps "
                "propellant name at a time."
            )

        if not np.isclose(composition_values[0], 1.0, rtol=0.0, atol=1e-6):
            raise ValueError(
                f"{self.name}: single propellant fraction must be 1.0. "
                f"Got {composition_values[0]}."
            )

        new_propellant = self._propellant_argument_from_composition()

        # Species changed, so rebuild the RocketProps backend.
        if new_propellant != self._rocketprops_propellant:
            self._rocketprops_propellant = new_propellant

            self._Propellant = Propellant(
                self._rocketprops_propellant,
                **{
                    name: getattr(self, name).value
                    for name in self._flash_names
                },
            )

            self._last_composition_values = composition_values
            self._last_flash_values = None
            self._property_cache.clear()
            return

        self._last_composition_values = composition_values
        self._last_flash_values = None
        self._property_cache.clear()

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

    @property
    def ignored_export_attributes(self) -> set[str]:
        return super().ignored_export_attributes | {
            "property_states",
            "_property_states",
            "external_property_names",
            "_external_property_names",
            "flash_names",
            "_flash_names",
            "Propellant",
            "_Propellant",
            "input_map",
            "_input_map",
            "rocketprops_propellant",
            "_rocketprops_propellant",
            "last_flash_values",
            "_last_flash_values",
            "property_cache",
            "_property_cache",
            "composition",
            "last_composition_values",
            "_last_composition_values",
        }