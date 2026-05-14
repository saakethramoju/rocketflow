from __future__ import annotations

from typing import Dict, List, Tuple, Union

import numpy as np
from scipy.optimize import root_scalar
import pyromat as pm


class IdealGas:
    """
    PYroMat ideal-gas wrapper with a Fluid-like API.

    Units:
        P   Pa
        T   K
        h   J/kg
        rho kg/m^3
        cp  J/kg-K
        cv  J/kg-K
        s   J/kg-K
        R   J/kg-K

    Supports:
        IdealGas("Nitrogen", T=...)
        IdealGas("Nitrogen", h=...)
        IdealGas("Nitrogen", P=..., T=...)
        IdealGas("Nitrogen", P=..., h=...)
        IdealGas({"Nitrogen": 0.78, "Oxygen": 0.21}, basis="mole", T=...)
        IdealGas({"Nitrogen": 0.75, "Oxygen": 0.25}, basis="mass", h=...)
    """

    _ALIASES = {
        "air": "ig.air",
        "n2": "ig.N2",
        "gn2": "ig.N2",
        "nitrogen": "ig.N2",
        "o2": "ig.O2",
        "oxygen": "ig.O2",
        "go2": "ig.O2",
        "h2": "ig.H2",
        "hydrogen": "ig.H2",
        "he": "ig.He",
        "helium": "ig.He",
        "ar": "ig.Ar",
        "argon": "ig.Ar",
        "co2": "ig.CO2",
        "carbon-dioxide": "ig.CO2",
        "carbon dioxide": "ig.CO2",
        "co": "ig.CO",
        "carbon-monoxide": "ig.CO",
        "carbon monoxide": "ig.CO",
        "ch4": "ig.CH4",
        "methane": "ig.CH4",
        "h2o": "ig.H2O",
        "water": "ig.H2O",
        "steam": "ig.H2O",
    }

    _RU = 8.31446261815324  # J/mol-K

    def __init__(
        self,
        fluid: Union[str, Dict[str, float]],
        basis: str = "mass",
        P: float | None = None,
        h: float | None = None,
        T: float | None = None,
        Q: float | None = None,
    ):
        self._configure_units()

        if Q is not None:
            raise ValueError("IdealGas does not support vapor quality Q.")

        self._species_ids: List[str] = []
        self._display_names: List[str] = []

        if isinstance(fluid, str):
            sid, display = self._normalize_name(fluid)
            self._species_ids = [sid]
            self._display_names = [display]
            self._mole_fractions = np.array([1.0])
            self._mass_fractions = np.array([1.0])
            self._mixture = False

        elif isinstance(fluid, dict):
            if basis not in ("mass", "mole"):
                raise ValueError("basis must be 'mass' or 'mole'")

            tmp: Dict[str, Tuple[float, List[str]]] = {}
            for user_name, frac in fluid.items():
                sid, display = self._normalize_name(user_name)
                total, names = tmp.get(sid, (0.0, []))
                tmp[sid] = (total + float(frac), names + [display])

            self._species_ids = list(tmp.keys())
            fractions = np.array([v[0] for v in tmp.values()], dtype=float)
            self._display_names = [", ".join(sorted(set(v[1]))) for v in tmp.values()]

            if not np.isclose(fractions.sum(), 1.0, atol=1e-6):
                raise ValueError(f"{basis.capitalize()} fractions must sum to 1.0")

            if basis == "mole":
                self._mole_fractions = fractions
                self._mass_fractions = self.mole_to_mass(self._species_ids, fractions)
            else:
                self._mass_fractions = fractions
                self._mole_fractions = self.mass_to_mole(self._species_ids, fractions)

            self._mixture = len(self._species_ids) > 1

        else:
            raise TypeError("fluid must be a string or a dict mixture")

        self._species = [pm.get(sid) for sid in self._species_ids]

        self._P: float | None = None
        self._h: float | None = None
        self._T: float | None = None

        if T is not None:
            self._P = None if P is None else float(P)
            self._T = float(T)
            self._h = self._enthalpy_from_T(self._T)

        elif h is not None:
            self._P = None if P is None else float(P)
            self._h = float(h)
            self._T = self._temperature_from_h(self._h)

        else:
            raise LookupError("Please provide either T or h.")

    # ---------------- Units ---------------- #

    @staticmethod
    def _configure_units():
        pm.config["unit_pressure"] = "Pa"
        pm.config["unit_temperature"] = "K"
        pm.config["unit_energy"] = "J"
        pm.config["unit_matter"] = "kg"
        pm.config["unit_volume"] = "m3"
        pm.config["unit_molar"] = "mol"

    # ---------------- Internal helpers ---------------- #

    def _require_pressure(self, property_name: str = "This property"):
        if self._P is None:
            raise ValueError(f"{property_name} requires pressure. Set gas.pressure first.")

    def _mix_mass_weighted(self, method: str, *, T: float | None = None, p: float | None = None):
        vals = []
        for sp in self._species:
            fn = getattr(sp, method)
            kwargs = {}
            if T is not None:
                kwargs["T"] = T
            if p is not None:
                kwargs["p"] = p
            vals.append(float(np.asarray(fn(**kwargs)).squeeze()))
        return float(np.dot(self._mass_fractions, vals))

    def _enthalpy_from_T(self, T: float) -> float:
        return self._mix_mass_weighted("h", T=T)

    def _temperature_from_h(self, h_target: float) -> float:
        def residual(T):
            return self._enthalpy_from_T(T) - h_target

        Tmin = self.minimum_temperature
        Tmax = self.maximum_temperature

        Ts = np.linspace(Tmin, Tmax, 400)
        rs = np.array([residual(T) for T in Ts])

        for T, r in zip(Ts, rs):
            if abs(r) < 1e-8:
                return float(T)

        for T1, T2, r1, r2 in zip(Ts[:-1], Ts[1:], rs[:-1], rs[1:]):
            if np.isfinite(r1) and np.isfinite(r2) and r1 * r2 <= 0:
                sol = root_scalar(residual, bracket=(T1, T2), method="brentq")
                return float(sol.root)

        raise ValueError(
            f"Could not solve ideal-gas T from h={h_target:.6g} J/kg "
            f"over T=[{Tmin:.3f}, {Tmax:.3f}] K."
        )

    def _partial_pressures(self) -> np.ndarray:
        self._require_pressure("Partial pressures")
        return self._mole_fractions * self._P

    # ---------------- Fractions ---------------- #

    @property
    def mole_fractions(self) -> dict:
        return {
            name: float(x)
            for name, x in zip(self._display_names, self._mole_fractions)
        }

    @mole_fractions.setter
    def mole_fractions(self, value: List[float]):
        if len(self._species_ids) == 1:
            raise ValueError("Cannot change mole fractions for a pure gas")
        if not np.isclose(sum(value), 1.0, atol=1e-6):
            raise ValueError("Mole fractions must sum to 1.0")
        self._mole_fractions = np.array(value, dtype=float)
        self._mass_fractions = self.mole_to_mass(self._species_ids, value)
        if self._h is not None:
            self._T = self._temperature_from_h(self._h)

    @property
    def mass_fractions(self) -> dict:
        return {
            name: float(x)
            for name, x in zip(self._display_names, self._mass_fractions)
        }

    @mass_fractions.setter
    def mass_fractions(self, value: List[float]):
        if len(self._species_ids) == 1:
            raise ValueError("Cannot change mass fractions for a pure gas")
        if not np.isclose(sum(value), 1.0, atol=1e-6):
            raise ValueError("Mass fractions must sum to 1.0")
        self._mass_fractions = np.array(value, dtype=float)
        self._mole_fractions = self.mass_to_mole(self._species_ids, value)
        if self._h is not None:
            self._T = self._temperature_from_h(self._h)

    # ---------------- State setters ---------------- #

    @property
    def pressure(self) -> float | None:
        return self._P

    @pressure.setter
    def pressure(self, value: float):
        self._P = float(value)

    @property
    def enthalpy(self) -> float:
        return self._h

    @enthalpy.setter
    def enthalpy(self, value: float):
        self._h = float(value)
        self._T = self._temperature_from_h(self._h)

    @property
    def temperature(self) -> float:
        return self._T

    @temperature.setter
    def temperature(self, value: float):
        self._T = float(value)
        self._h = self._enthalpy_from_T(self._T)

    @property
    def HP(self) -> Tuple[float, float | None]:
        return self._h, self._P

    @HP.setter
    def HP(self, values: Tuple[float, float]):
        if not isinstance(values, (tuple, list)) or len(values) != 2:
            raise ValueError("HP must be set with (h, P)")
        self._h = float(values[0])
        self._P = None if values[1] is None else float(values[1])
        self._T = self._temperature_from_h(self._h)

    @property
    def TP(self) -> Tuple[float, float | None]:
        return self._T, self._P

    @TP.setter
    def TP(self, values: Tuple[float, float]):
        if not isinstance(values, (tuple, list)) or len(values) != 2:
            raise ValueError("TP must be set with (T, P)")
        self._T = float(values[0])
        self._P = None if values[1] is None else float(values[1])
        self._h = self._enthalpy_from_T(self._T)

    # ---------------- Thermo properties ---------------- #

    @property
    def species(self) -> List[str]:
        return self._display_names

    @property
    def phase(self) -> str:
        return "Ideal Gas"

    @property
    def compressibility(self) -> float:
        return 1.0

    @property
    def density(self) -> float:
        self._require_pressure("Density")
        return self._P / (self.gas_constant * self._T)

    @property
    def specific_volume(self) -> float:
        return 1.0 / self.density

    @property
    def molar_mass(self) -> float:
        return float(1.0 / np.sum(self._mass_fractions / self._molar_masses()))

    @property
    def gas_constant(self) -> float:
        return self._RU / self.molar_mass

    @property
    def specific_heat_cp(self) -> float:
        return self._mix_mass_weighted("cp", T=self._T)

    @property
    def specific_heat_cv(self) -> float:
        return self._mix_mass_weighted("cv", T=self._T)

    @property
    def specific_heat(self) -> float:
        return self.specific_heat_cp

    @property
    def specific_heat_ratio(self) -> float:
        cv = self.specific_heat_cv
        return None if cv == 0 else self.specific_heat_cp / cv

    @property
    def internal_energy(self) -> float:
        return self._mix_mass_weighted("e", T=self._T)

    @property
    def free_energy(self) -> float:
        self._require_pressure("Free energy")
        try:
            return self._mix_mass_weighted("f", T=self._T, p=self._P)
        except Exception:
            return self.internal_energy - self._T * self.entropy

    @property
    def gibbs_energy(self) -> float:
        self._require_pressure("Gibbs energy")
        try:
            if not self._mixture:
                return self._mix_mass_weighted("g", T=self._T, p=self._P)

            vals = []
            for wi, sp, pi in zip(self._mass_fractions, self._species, self._partial_pressures()):
                vals.append(wi * float(np.asarray(sp.g(T=self._T, p=pi)).squeeze()))
            return float(sum(vals))

        except Exception:
            return self.enthalpy - self._T * self.entropy

    @property
    def entropy(self) -> float:
        self._require_pressure("Entropy")

        if not self._mixture:
            return self._mix_mass_weighted("s", T=self._T, p=self._P)

        vals = []
        for wi, sp, pi in zip(self._mass_fractions, self._species, self._partial_pressures()):
            vals.append(wi * float(np.asarray(sp.s(T=self._T, p=pi)).squeeze()))
        return float(sum(vals))

    @property
    def quality(self) -> float:
        return 1.0

    @quality.setter
    def quality(self, value: float):
        raise ValueError("IdealGas does not support vapor quality.")

    @property
    def speed_of_sound(self) -> float:
        return float(np.sqrt(self.specific_heat_ratio * self.gas_constant * self._T))

    @property
    def minimum_pressure(self) -> float:
        return 1e-9

    @property
    def maximum_pressure(self) -> float:
        return np.inf

    @property
    def minimum_temperature(self) -> float:
        mins = []
        for sp in self._species:
            try:
                mins.append(float(sp.Tlim()[0]))
            except Exception:
                mins.append(200.0)
        return max(mins)

    @property
    def maximum_temperature(self) -> float:
        maxs = []
        for sp in self._species:
            try:
                maxs.append(float(sp.Tlim()[1]))
            except Exception:
                maxs.append(6000.0)
        return min(maxs)

    @property
    def is_mixture(self) -> bool:
        return self._mixture

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
            ("Gas(es)", ", ".join(self._display_names)),
            ("Mole fractions", format_dict(self.mole_fractions, 3)),
            ("Mass fractions", format_dict(self.mass_fractions, 3)),
            ("Phase", self.phase),
            ("Pressure [Pa]", self._safe(self.pressure, ".3e")),
            ("Temperature [K]", self._safe(self.temperature, ".2f")),
            ("Density [kg/m³]", self._safe(self.density, ".3f") if self._P is not None else "N/A"),
            ("Compressibility Z", self._safe(self.compressibility, ".3f")),
            ("Internal energy [J/kg]", self._safe(self.internal_energy, ".3e")),
            ("Enthalpy [J/kg]", self._safe(self.enthalpy, ".3e")),
            ("Entropy [J/kg-K]", self._safe(self.entropy, ".3e") if self._P is not None else "N/A"),
            ("Cp [J/kg-K]", self._safe(self.specific_heat_cp, ".3f")),
            ("Cv [J/kg-K]", self._safe(self.specific_heat_cv, ".3f")),
            ("Gamma", self._safe(self.specific_heat_ratio, ".5f")),
            ("Gas constant [J/kg-K]", self._safe(self.gas_constant, ".3f")),
            ("Molar mass [kg/mol]", self._safe(self.molar_mass, ".6f")),
            ("Speed of sound [m/s]", self._safe(self.speed_of_sound, ".3f")),
        ]

        width = max(len(r[0]) for r in rows)
        return "\n".join(f"{key:<{width}} : {val}" for key, val in rows)

    def __repr__(self) -> str:
        species_str = ", ".join(self._display_names)
        P_str = "None" if self._P is None else f"{self._P:.3e}"
        return (
            f"{self.__class__.__name__}(species=[{species_str}], "
            f"P={P_str} Pa, h={self._h:.3e} J/kg, T={self.temperature:.2f} K)"
        )

    # ---------------- Utilities ---------------- #

    @staticmethod
    def _alias_key(name: str) -> str:
        return name.strip().lower().replace("_", "-")

    @classmethod
    def _normalize_name(cls, user_name: str) -> Tuple[str, str]:

        key = cls._alias_key(user_name)

        sid = cls._ALIASES.get(key, user_name)

        if not sid.startswith("ig."):
            sid = f"ig.{sid}"

        try:
            pm.get(sid)
        except Exception:
            raise ValueError(
                f"Invalid ideal gas '{user_name}'. "
                f"Use IdealGas.show_available_gases() to check valid names."
            )

        return sid, user_name

    @classmethod
    def add_alias(cls, alias: str, pyromat_name: str) -> None:
        if not pyromat_name.startswith("ig."):
            pyromat_name = f"ig.{pyromat_name}"
        cls._ALIASES[cls._alias_key(alias)] = pyromat_name

    @classmethod
    def add_aliases(cls, aliases: dict[str, str]) -> None:
        for alias, pyromat_name in aliases.items():
            cls.add_alias(alias, pyromat_name)

    @classmethod
    def remove_alias(cls, alias: str) -> None:
        cls._ALIASES.pop(cls._alias_key(alias), None)

    @classmethod
    def show_aliases(cls) -> dict[str, str]:
        width = max(len(alias) for alias in cls._ALIASES)
        print("IdealGas Aliases")
        print("-" * (width + 20))
        for alias, backend in sorted(cls._ALIASES.items()):
            print(f"{alias:<{width}} -> {backend}")
        return dict(cls._ALIASES)

    @staticmethod
    def get_available_gases() -> List[str]:
        gases = sorted(k for k in pm.dat.data.keys() if k.startswith("ig."))
        return [g[3:] for g in gases]

    @staticmethod
    def show_available_gases() -> List[str]:
        gases = IdealGas.get_available_gases()
        for g in gases:
            print(g)
        return gases

    @staticmethod
    def _molar_mass_of(species_id: str) -> float:
        if not species_id.startswith("ig."):
            species_id = f"ig.{species_id}"
        return float(np.asarray(pm.get(species_id).mw()).squeeze())

    def _molar_masses(self) -> np.ndarray:
        return np.array([self._molar_mass_of(sid) for sid in self._species_ids], dtype=float)

    @staticmethod
    def mole_to_mass(species_ids: List[str], mole_fractions: List[float]):
        if not np.isclose(sum(mole_fractions), 1.0, atol=1e-6):
            raise ValueError("Mole fractions must sum to 1.0")
        x = np.asarray(mole_fractions, dtype=float)
        M = np.array([IdealGas._molar_mass_of(sid) for sid in species_ids])
        return x * M / np.dot(x, M)

    @staticmethod
    def mass_to_mole(species_ids: List[str], mass_fractions: List[float]):
        if not np.isclose(sum(mass_fractions), 1.0, atol=1e-6):
            raise ValueError("Mass fractions must sum to 1.0")
        w = np.asarray(mass_fractions, dtype=float)
        M = np.array([IdealGas._molar_mass_of(sid) for sid in species_ids])
        inv = w / M
        return inv / inv.sum()


if __name__ == "__main__":
    gas = IdealGas("Nitrogen", T=300)
    print(gas)

    print("-------------------")

    gas.pressure = 101325
    print(gas)

    print("-------------------")

    air = IdealGas({"Nitrogen": 0.78, "Oxygen": 0.21, "Argon": 0.01}, basis="mole", T=298.15)
    print(air)

    print("-------------------")

    same_air = IdealGas({"Nitrogen": 0.78, "Oxygen": 0.21, "Argon": 0.01}, basis="mole", h=air.enthalpy)
    print(same_air)