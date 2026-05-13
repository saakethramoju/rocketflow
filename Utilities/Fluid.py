from typing import List, Union, Dict, Tuple

import numpy as np
from scipy.optimize import root_scalar

import CoolProp.CoolProp as CP


class Fluid:
    """
    High-level wrapper around CoolProp AbstractState objects with added
    convenience methods for initialization, state updates, and property access.

    Supports both pure fluids and mixtures by mole or mass fractions.

    Notes
    -----
    Public API is intentionally kept the same as the old pyfluids-backed class:
    Fluid(..., basis="mass", P=None, h=None, T=None, Q=None), HP, pressure,
    enthalpy, temperature, density, quality, etc.
    """

    _ALIASES = {
        "rp-1": "n-Dodecane",
        "rp1": "n-Dodecane",
        "jet-a": "n-Dodecane",
        "jeta": "n-Dodecane",
        "kerosene": "n-Dodecane",
        "lox": "Oxygen",
        "water": "Water",
        "air": "Air",
        "nitrous": "NitrousOxide",
        "n2o": "NitrousOxide",
        "gn2": "Nitrogen",
        "n2": "Nitrogen"
    }

    _PHASE_NAMES = {
        getattr(CP, "iphase_unknown", -999): "Unknown",
        getattr(CP, "iphase_liquid", -999): "Liquid",
        getattr(CP, "iphase_supercritical", -999): "Supercritical",
        getattr(CP, "iphase_supercritical_gas", -999): "SupercriticalGas",
        getattr(CP, "iphase_supercritical_liquid", -999): "SupercriticalLiquid",
        getattr(CP, "iphase_gas", -999): "Gas",
        getattr(CP, "iphase_twophase", -999): "TwoPhase",
        getattr(CP, "iphase_critical_point", -999): "CriticalPoint",
    }

    def __init__(
        self,
        fluid: Union[str, Dict[str, float]],
        basis: str = "mass",
        P: float = None,
        h: float = None,
        T: float = None,
        Q: float = None,
    ):
        """
        Initialize a Fluid state.

        Parameters
        ----------
        fluid : str or dict
            Pure fluid name or dictionary of {fluid_name: fraction}.
        basis : str
            "mole" or "mass" basis for mixture fractions. Default is "mass".
        P : float, optional
            Absolute pressure in Pa.
        h : float, optional
            Specific enthalpy in J/kg.
        T : float, optional
            Absolute temperature in K.
        Q : float, optional
            Vapor quality from 0 to 1.
        """
        valid_fluids = Fluid.get_available_fluids()

        self._fluids: List[str] = []
        self._display_names: List[str] = []

        if isinstance(fluid, str):
            backend, display = Fluid._normalize_name(fluid)
            if backend not in valid_fluids:
                raise ValueError(
                    f"Invalid fluid '{fluid}'. "
                    f"Use Fluid.show_available_fluids() to check valid names."
                )
            self._fluids = [backend]
            self._display_names = [display]
            self._mole_fractions = np.array([1.0])
            self._mass_fractions = np.array([1.0])
            self._mixture = False

        elif isinstance(fluid, dict):
            if len(fluid) == 1:
                f, frac = next(iter(fluid.items()))
                if not np.isclose(frac, 1.0, atol=1e-12):
                    raise ValueError(f"Single-component dict must have fraction = 1.0, got {frac}")
                backend, display = Fluid._normalize_name(f)
                if backend not in valid_fluids:
                    raise ValueError(
                        f"Invalid fluid '{f}'. "
                        f"Use Fluid.show_available_fluids() to check valid names."
                    )
                self._fluids = [backend]
                self._display_names = [display]
                self._mole_fractions = np.array([1.0])
                self._mass_fractions = np.array([1.0])
                self._mixture = False

            else:
                tmp: Dict[str, Tuple[float, List[str]]] = {}
                for user_name, frac in fluid.items():
                    backend, display = Fluid._normalize_name(user_name)
                    if backend not in valid_fluids:
                        raise ValueError(
                            f"Invalid fluid '{user_name}' (backend '{backend}' not found). "
                            f"Use Fluid.show_available_fluids() to check valid names."
                        )
                    total, names = tmp.get(backend, (0.0, []))
                    tmp[backend] = (total + float(frac), names + [display])

                self._fluids = list(tmp.keys())
                fractions = np.array([v[0] for v in tmp.values()], dtype=float)
                self._display_names = [", ".join(sorted(set(names))) for _, names in tmp.values()]

                if basis == "mole":
                    if not np.isclose(fractions.sum(), 1.0, atol=1e-6):
                        raise ValueError("Mole fractions must sum to 1.0")
                    self._mole_fractions = fractions
                    self._mass_fractions = Fluid.mole_to_mass(self._fluids, fractions)
                elif basis == "mass":
                    if not np.isclose(fractions.sum(), 1.0, atol=1e-6):
                        raise ValueError("Mass fractions must sum to 1.0")
                    self._mass_fractions = fractions
                    self._mole_fractions = Fluid.mass_to_mole(self._fluids, fractions)
                else:
                    raise ValueError("basis must be 'mole' or 'mass'")

                self._mixture = len(self._fluids) > 1

        else:
            raise TypeError("fluid must be a string (pure) or dict (mixture)")

        self._P = None
        self._h = None
        self._fluid_string = "&".join(self._fluids)
        self._backend = self._build_state()
        self._pyfluid = self._backend  # backward-compatible internal alias

        if P is not None and h is not None:
            self._P, self._h = float(P), float(h)
        elif P is not None and T is not None:
            self._P = float(P)
            self._h = self._enthalpy_from_PT(float(P), float(T))
        elif P is not None and Q is not None:
            self._P = float(P)
            self._h = self._enthalpy_from_PQ(float(P), float(Q))
        elif T is not None and Q is not None:
            self._P, self._h = self._state_from_TQ(float(T), float(Q))
        else:
            raise LookupError("Please provide at least two thermodynamic properties!")

        self.set_pyfluid()

    # ---------------- Core ---------------- #
    def _build_state(self):
        """Create and configure a CoolProp AbstractState."""
        state = CP.AbstractState("HEOS", self._fluid_string)
        if self._mixture:
            state.set_mass_fractions([float(x) for x in self._mass_fractions])
        return state

    def _update_state(self, input_pair, value1: float, value2: float):
        """Update CoolProp state and keep the old internal alias current."""
        self._backend.update(input_pair, float(value1), float(value2))
        self._pyfluid = self._backend

    def set_pyfluid(self):
        """Rebuild backend CoolProp state using current P and h."""
        if self._mixture:
            T, Q = Fluid.get_temperature_and_quality(self._backend, self._P, self._h)
            if 0.0 < Q < 1.0:
                self._update_state(CP.PQ_INPUTS, self._P, Q)
            else:
                self._update_state(CP.PT_INPUTS, self._P, T)
        else:
            self._update_state(CP.HmassP_INPUTS, self._h, self._P)

    # ---------------- Internal state helpers ---------------- #
    def _enthalpy_from_PT(self, P: float, T: float) -> float:
        self._update_state(CP.PT_INPUTS, P, T)
        return float(self._backend.hmass())

    def _enthalpy_from_PQ(self, P: float, Q: float) -> float:
        self._update_state(CP.PQ_INPUTS, P, Q)
        return float(self._backend.hmass())

    def _state_from_TQ(self, T: float, Q: float) -> Tuple[float, float]:
        self._update_state(CP.QT_INPUTS, Q, T)
        return float(self._backend.p()), float(self._backend.hmass())

    def _keyed_output(self, key, default=None):
        try:
            return float(self._backend.keyed_output(key))
        except Exception:
            return default

    def _trivial_output(self, key, default=None):
        try:
            return float(self._backend.trivial_keyed_output(key))
        except Exception:
            return default

    # ---------------- Fractions ---------------- #
    @property
    def mole_fractions(self) -> dict:
        """Return mole fractions as {fluid_name: value}."""
        return {f: float(x) for f, x in zip(self._fluids, self._mole_fractions)}

    @mole_fractions.setter
    def mole_fractions(self, value: List[float]):
        """Update mole fractions. Fractions must sum to 1."""
        if len(self._fluids) == 1:
            raise ValueError("Cannot change mole fractions for a pure fluid")
        if not np.isclose(sum(value), 1.0, atol=1e-6):
            raise ValueError("Mole fractions must sum to 1.0")
        self._mole_fractions = np.array(value, dtype=float)
        self._mass_fractions = Fluid.mole_to_mass(self._fluids, value)
        self._backend = self._build_state()
        self._pyfluid = self._backend
        if self._P is not None and self._h is not None:
            self.set_pyfluid()

    @property
    def mass_fractions(self) -> dict:
        """Return mass fractions as {fluid_name: value}."""
        return {f: float(x) for f, x in zip(self._fluids, self._mass_fractions)}

    @mass_fractions.setter
    def mass_fractions(self, value: List[float]):
        """Update mass fractions. Fractions must sum to 1."""
        if len(self._fluids) == 1:
            raise ValueError("Cannot change mass fractions for a pure fluid")
        if not np.isclose(sum(value), 1.0, atol=1e-6):
            raise ValueError("Mass fractions must sum to 1.0")
        self._mass_fractions = np.array(value, dtype=float)
        self._mole_fractions = Fluid.mass_to_mole(self._fluids, value)
        self._backend = self._build_state()
        self._pyfluid = self._backend
        if self._P is not None and self._h is not None:
            self.set_pyfluid()

    # ---------------- State setters ---------------- #
    @property
    def pressure(self) -> float:
        """Absolute pressure in Pa."""
        return self._P

    @pressure.setter
    def pressure(self, value: float):
        self._P = float(value)
        if self._h is not None:
            self.set_pyfluid()

    @property
    def enthalpy(self) -> float:
        """Mass-specific enthalpy in J/kg."""
        return self._h

    @enthalpy.setter
    def enthalpy(self, value: float):
        self._h = float(value)
        if self._P is not None:
            self.set_pyfluid()

    @property
    def HP(self) -> Tuple[float, float]:
        """Return (h [J/kg], P [Pa])."""
        return self._h, self._P

    @HP.setter
    def HP(self, values: Tuple[float, float]):
        """Update enthalpy and pressure simultaneously."""
        if not isinstance(values, (tuple, list)) or len(values) != 2:
            raise ValueError("HP must be set with (P, h)")
        self._h, self._P = float(values[0]), float(values[1])
        self.set_pyfluid()

    @property
    def TP(self) -> Tuple[float, float]:
        """Return (T [K], P [Pa])."""
        return self.temperature, self._P


    @TP.setter
    def TP(self, values: Tuple[float, float]):
        """Update temperature and pressure simultaneously."""
        if not isinstance(values, (tuple, list)) or len(values) != 2:
            raise ValueError("TP must be set with (T, P)")

        T = float(values[0])
        P = float(values[1])

        self._P = P
        self._h = self._enthalpy_from_PT(P, T)

        self.set_pyfluid()

    # ---------------- Thermo properties ---------------- #
    @property
    def species(self) -> List[str]:
        return self._display_names

    @property
    def temperature(self) -> float:
        """Absolute temperature in K."""
        return float(self._backend.T())
    
    @temperature.setter
    def temperature(self, value: float):
        """
        Update temperature while holding pressure constant.

        Requires pressure to already be defined.
        """
        if self._P is None:
            raise ValueError(
                "Cannot set temperature without pressure. "
                "Set pressure first."
            )

        T = float(value)

        self._h = self._enthalpy_from_PT(self._P, T)
        self.set_pyfluid()

    @property
    def phase(self) -> str:
        """Thermodynamic phase name, mapped to the old pyfluids-style names."""
        try:
            return Fluid._PHASE_NAMES.get(int(self._backend.phase()), "Unknown")
        except Exception:
            return "Unknown"

    @property
    def compressibility(self) -> float:
        """Compressibility factor Z."""
        return self._keyed_output(CP.iZ)

    @property
    def conductivity(self) -> float:
        """Thermal conductivity in W/m-K."""
        try:
            return float(self._backend.conductivity())
        except Exception:
            return None

    @property
    def critical_pressure(self) -> float:
        """Critical pressure in Pa."""
        try:
            return float(self._backend.p_critical())
        except Exception:
            return self._trivial_output(CP.iP_critical)

    @property
    def critical_temperature(self) -> float:
        """Critical temperature in K."""
        try:
            return float(self._backend.T_critical())
        except Exception:
            return self._trivial_output(CP.iT_critical)

    @property
    def density(self) -> float:
        """Mass density in kg/m^3."""
        return float(self._backend.rhomass())

    @property
    def dynamic_viscosity(self) -> float:
        """Dynamic viscosity in Pa-s."""
        try:
            return float(self._backend.viscosity())
        except Exception:
            return None

    @property
    def entropy(self) -> float:
        """Mass-specific entropy in J/kg-K."""
        return float(self._backend.smass())

    @property
    def freezing_temperature(self) -> float:
        """Freezing/melting temperature in K when available; otherwise Tmin."""
        # CoolProp does not expose a universal pyfluids-equivalent freezing point.
        # Tmin is the most useful direct replacement for bounds checking.
        return self.minimum_temperature

    @property
    def internal_energy(self) -> float:
        """Mass-specific internal energy in J/kg."""
        return float(self._backend.umass())

    @property
    def kinematic_viscosity(self) -> float:
        """Kinematic viscosity in m^2/s."""
        mu = self.dynamic_viscosity
        rho = self.density
        if mu is None or rho is None or rho == 0:
            return None
        return mu / rho

    @property
    def maximum_pressure(self) -> float:
        """Maximum valid pressure in Pa."""
        try:
            return float(self._backend.pmax())
        except Exception:
            return self._trivial_output(CP.iP_max)

    @property
    def maximum_temperature(self) -> float:
        """Maximum valid temperature in K."""
        try:
            return float(self._backend.Tmax())
        except Exception:
            return self._trivial_output(CP.iT_max)

    @property
    def minimum_pressure(self) -> float:
        """Minimum valid pressure in Pa."""
        try:
            return float(self._backend.p_triple())
        except Exception:
            return self._trivial_output(CP.iP_triple)

    @property
    def minimum_temperature(self) -> float:
        """Minimum valid temperature in K."""
        try:
            return float(self._backend.Tmin())
        except Exception:
            return self._trivial_output(CP.iT_min)

    @property
    def molar_mass(self) -> float:
        """Molar mass in kg/mol."""
        return float(self._backend.molar_mass())

    @property
    def prandtl(self) -> float:
        """Prandtl number."""
        return self._keyed_output(CP.iPrandtl)

    @property
    def speed_of_sound(self) -> float:
        """Speed of sound in m/s."""
        try:
            return float(self._backend.speed_sound())
        except Exception:
            return None
    @property
    def specific_heat(self) -> float:
        """Mass-specific heat capacity at constant pressure, Cp, in J/kg-K."""
        try:
            return float(self._backend.cpmass())
        except Exception:
            return None

    @property
    def specific_volume(self) -> float:
        """Specific volume in m^3/kg."""
        rho = self.density
        if rho is None or rho == 0:
            return None
        return 1.0 / rho

    @property
    def surface_tension(self) -> float:
        """Surface tension in N/m when available."""
        try:
            return float(self._backend.surface_tension())
        except Exception:
            return None

    @property
    def triple_pressure(self) -> float:
        """Triple point pressure in Pa."""
        try:
            return float(self._backend.p_triple())
        except Exception:
            return self._trivial_output(CP.iP_triple)

    @property
    def triple_temperature(self) -> float:
        """Triple point temperature in K."""
        try:
            return float(self._backend.Ttriple())
        except Exception:
            return self._trivial_output(CP.iT_triple)

    @property
    def is_mixture(self) -> bool:
        """Return True if this fluid is a mixture, False if pure."""
        return self._mixture

    @property
    def quality(self) -> float:
        """
        Vapor quality from 0 to 1. Only physically meaningful in TwoPhase.

        Returns 0.0 for liquid-like single phase and 1.0 for gas/supercritical
        single phase, matching the convention in the old class.
        """
        ph = self.phase
        if ph == "TwoPhase":
            try:
                return float(self._backend.Q())
            except Exception:
                return float("nan")
        if ph in ("Gas", "Supercritical", "SupercriticalGas"):
            return 1.0
        if ph in ("Liquid", "SupercriticalLiquid"):
            return 0.0
        return float("nan")
    
    @quality.setter
    def quality(self, value: float):
        """
        Update vapor quality while holding pressure constant.

        Requires pressure to already be defined.
        """
        if self._P is None:
            raise ValueError(
                "Cannot set quality without pressure. "
                "Set pressure first."
            )

        Q = float(value)

        if not (0.0 <= Q <= 1.0):
            raise ValueError("Quality must be between 0 and 1.")

        self._h = self._enthalpy_from_PQ(self._P, Q)
        self.set_pyfluid()

    @property
    def saturation_temperature(self) -> float:
        """Saturation temperature in K for current pressure, only if P <= Pc."""
        pc = self.critical_pressure
        if pc is not None and self.pressure > pc:
            return None
        try:
            tmp = self._build_state()
            tmp.update(CP.PQ_INPUTS, self._P, 1.0)
            return float(tmp.T())
        except Exception:
            return None

    # ---------------- String output ---------------- #
    def _safe(self, value, fmt=".3e"):
        if value is None:
            return "N/A"
        try:
            return f"{value:{fmt}}"
        except Exception:
            return str(value)

    def __str__(self):
        def format_dict(d: dict, decimals=3):
            return {k: round(v, decimals) for k, v in d.items()}

        rows = [
            ("Fluid(s)", ", ".join(self._display_names)),
            ("Mole fractions", format_dict(self.mole_fractions, 3)),
            ("Mass fractions", format_dict(self.mass_fractions, 3)),
            ("Phase", self.phase),
            ("Pressure [Pa]", self._safe(self.pressure, ".3e")),
            ("Temperature [K]", self._safe(self.temperature, ".2f")),
            ("Density [kg/m³]", self._safe(self.density, ".3f")),
            ("Quality", self._safe(self.quality, ".3f")),
            ("Internal energy [J/kg]", self._safe(self.internal_energy, ".3e")),
            ("Enthalpy [J/kg]", self._safe(self.enthalpy, ".3e")),
            ("Entropy [J/kg-K]", self._safe(self.entropy, ".3e")),
            ("Dynamic viscosity [Pa·s]", self._safe(self.dynamic_viscosity, ".3e")),
            ("Conductivity [W/m-K]", self._safe(self.conductivity, ".3f")),
            ("Saturation temperature [K]", self._safe(self.saturation_temperature, ".2f")),
            ("Molar mass [kg/mol]", self._safe(self.molar_mass, ".6f")),
        ]
        width = max(len(r[0]) for r in rows)
        return "\n".join(f"{key:<{width}} : {val}" for key, val in rows)

    def __repr__(self) -> str:
        species_str = ", ".join(self._display_names)
        return (
            f"{self.__class__.__name__}(species=[{species_str}], "
            f"P={self._P:.3e} Pa, h={self._h:.3e} J/kg, T={self.temperature:.2f} K)"
        )

    # ---------------- Utilities ---------------- #
    @staticmethod
    def _alias_key(name: str) -> str:
        return name.strip().lower().replace(" ", "").replace("_", "-")


    @classmethod
    def _normalize_name(cls, user_name: str) -> Tuple[str, str]:
        """
        Return (backend_name, display_name). If user_name is an alias, map it to
        the CoolProp backend name but preserve user input for display.
        """
        display = user_name
        key = cls._alias_key(user_name)
        backend = cls._ALIASES.get(key, user_name)
        return backend, display


    @classmethod
    def add_alias(cls, alias: str, coolprop_name: str) -> None:
        cls._ALIASES[cls._alias_key(alias)] = coolprop_name


    @classmethod
    def add_aliases(cls, aliases: dict[str, str]) -> None:
        for alias, coolprop_name in aliases.items():
            cls.add_alias(alias, coolprop_name)


    @classmethod
    def remove_alias(cls, alias: str) -> None:
        cls._ALIASES.pop(cls._alias_key(alias), None)


    @classmethod
    def show_aliases(cls) -> dict[str, str]:

        width = max(len(alias) for alias in cls._ALIASES)

        print("Fluid Aliases")
        print("-" * (width + 20))

        for alias, backend in sorted(cls._ALIASES.items()):
            print(f"{alias:<{width}} -> {backend}")

        return dict(cls._ALIASES)

    @staticmethod
    def _molar_mass_of(fluid: str) -> float:
        """Return pure-fluid molar mass in kg/mol."""
        return float(CP.PropsSI("M", fluid))

    @staticmethod
    def mole_to_mass(fluids: List[str], mole_fractions: List[float]):
        """Convert mole fractions to mass fractions."""
        if not np.isclose(sum(mole_fractions), 1.0, atol=1e-6):
            raise ValueError("Mole fractions must sum to 1.0")
        mole_fractions = np.asarray(mole_fractions, dtype=float)
        molar_masses = np.array([Fluid._molar_mass_of(f) for f in fluids])
        m_bar = np.dot(mole_fractions, molar_masses)
        return mole_fractions * molar_masses / m_bar

    @staticmethod
    def mass_to_mole(fluids: List[str], mass_fractions: List[float]):
        """Convert mass fractions to mole fractions."""
        if not np.isclose(sum(mass_fractions), 1.0, atol=1e-6):
            raise ValueError("Mass fractions must sum to 1.0")
        mass_fractions = np.asarray(mass_fractions, dtype=float)
        molar_masses = np.array([Fluid._molar_mass_of(f) for f in fluids])
        inv = mass_fractions / molar_masses
        return inv / inv.sum()

    @staticmethod
    def get_temperature_and_quality(fluid, P: float, target_enthalpy: float) -> Tuple[float, float]:
        """
        Given a CoolProp AbstractState, pressure, and enthalpy, return (T, Q).

        For mixtures this avoids relying on direct P-H flashes through the dome.
        It mirrors the old pyfluids workaround: compare target enthalpy to the
        saturated liquid/vapor enthalpies, then solve T at fixed P outside dome.
        """
        try:
            fluid.update(CP.PQ_INPUTS, P, 0.0)
            h_liquid = float(fluid.hmass())
            T_sat = float(fluid.T())

            fluid.update(CP.PQ_INPUTS, P, 1.0)
            h_vapor = float(fluid.hmass())
        except Exception:
            h_liquid = None
            h_vapor = None
            T_sat = None

        h = float(target_enthalpy)

        if h_liquid is not None and h_vapor is not None and h_liquid <= h <= h_vapor:
            denom = h_vapor - h_liquid
            Q = 0.0 if abs(denom) < 1e-15 else (h - h_liquid) / denom
            return T_sat, float(Q)

        def residual(T):
            try:
                fluid.update(CP.PT_INPUTS, P, T)
                return float(fluid.hmass()) - h
            except Exception:
                return np.nan

        try:
            Tmin = float(fluid.Tmin())
        except Exception:
            Tmin = 1.0
        try:
            Tmax = float(fluid.Tmax())
        except Exception:
            Tmax = 5000.0

        # Build a robust bracket by scanning the valid temperature range.
        Ts = np.linspace(Tmin * 1.000001, Tmax * 0.999999, 300)
        vals = []
        for T in Ts:
            r = residual(T)
            vals.append(r if np.isfinite(r) else np.nan)

        bracket = None
        for T1, T2, r1, r2 in zip(Ts[:-1], Ts[1:], vals[:-1], vals[1:]):
            if not (np.isfinite(r1) and np.isfinite(r2)):
                continue
            if r1 == 0:
                bracket = (T1, T1)
                break
            if r1 * r2 <= 0:
                bracket = (T1, T2)
                break

        if bracket is None:
            raise ValueError(
                f"Could not find a valid temperature bracket for P={P:.6g} Pa, "
                f"h={h:.6g} J/kg over T=[{Tmin:.6g}, {Tmax:.6g}] K."
            )

        if bracket[0] == bracket[1]:
            T = bracket[0]
        else:
            sol = root_scalar(residual, method="brentq", bracket=bracket)
            T = float(sol.root)

        if h_liquid is not None and h < h_liquid:
            Q = 0.0
        elif h_vapor is not None and h > h_vapor:
            Q = 1.0
        else:
            Q = float("nan")

        return T, Q

    @staticmethod
    def show_available_fluids():
        """Print and return available CoolProp fluid names."""
        fluids = Fluid.get_available_fluids()
        for f in fluids:
            print(f)
        return fluids

    @staticmethod
    def get_available_fluids():
        """Return available CoolProp fluid names."""
        return sorted(CP.get_global_param_string("FluidsList").split(","))


if __name__ == "__main__":
    # f = Fluid({"Nitrogen": 0.78, "Oxygen": 0.21, "Argon": 0.01}, basis="mole", P=101325, T=298.15)
    # f = Fluid({"Nitrogen": 1}, P=101325, h=311200)
    # f = Fluid("Methane", P=3e6, Q=0.1)
    f = Fluid("LOX", P=450*6894, T=90)
    print(f)
    print(f.minimum_pressure)
    print(f.density)

    print("-------------------")
    f = Fluid("Air", P=101325, T=298.15)
    print(f)
