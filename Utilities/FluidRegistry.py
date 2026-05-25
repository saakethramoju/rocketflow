from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True)
class SpeciesRecord:
    name: str
    coolprop: str | None = None
    pyromat: str | None = None  # no "ig." prefix


SPECIES_DATABASE = MappingProxyType({
    "Air": SpeciesRecord("Air", coolprop="Air", pyromat="air"),
    "Argon": SpeciesRecord("Argon", coolprop="Argon", pyromat="Ar"),
    "CarbonDioxide": SpeciesRecord("CarbonDioxide", coolprop="CarbonDioxide", pyromat="CO2"),
    "CarbonMonoxide": SpeciesRecord("CarbonMonoxide", coolprop="CarbonMonoxide", pyromat="CO"),
    "Helium": SpeciesRecord("Helium", coolprop="Helium", pyromat="He"),
    "Hydrogen": SpeciesRecord("Hydrogen", coolprop="Hydrogen", pyromat="H2"),
    "Methane": SpeciesRecord("Methane", coolprop="Methane", pyromat="CH4"),
    "Nitrogen": SpeciesRecord("Nitrogen", coolprop="Nitrogen", pyromat="N2"),
    "Oxygen": SpeciesRecord("Oxygen", coolprop="Oxygen", pyromat="O2"),
    "Water": SpeciesRecord("Water", coolprop="Water", pyromat="H2O"),

    "Ammonia": SpeciesRecord("Ammonia", coolprop="Ammonia", pyromat="NH3"),
    "Ethane": SpeciesRecord("Ethane", coolprop="Ethane", pyromat="C2H6"),
    "Ethylene": SpeciesRecord("Ethylene", coolprop="Ethylene", pyromat="C2H4"),
    "Propane": SpeciesRecord("Propane", coolprop="n-Propane", pyromat="C3H8"),
    "Butane": SpeciesRecord("Butane", coolprop="Butane", pyromat="C4H10"),
    "IsoButane": SpeciesRecord("IsoButane", coolprop="IsoButane", pyromat="C4H10"),
    "Acetylene": SpeciesRecord("Acetylene", coolprop=None, pyromat="C2H2"),
    "Benzene": SpeciesRecord("Benzene", coolprop="Benzene", pyromat="C6H6"),
    "Toluene": SpeciesRecord("Toluene", coolprop="Toluene", pyromat="C7H8"),
    "Methanol": SpeciesRecord("Methanol", coolprop="Methanol", pyromat="CH4O"),
    "Ethanol": SpeciesRecord("Ethanol", coolprop="Ethanol", pyromat="C2H6O"),
    "NitrousOxide": SpeciesRecord("NitrousOxide", coolprop="NitrousOxide", pyromat="N2O"),
    "HydrogenChloride": SpeciesRecord("HydrogenChloride", coolprop="HydrogenChloride", pyromat="HCl"),
    "HydrogenSulfide": SpeciesRecord("HydrogenSulfide", coolprop="HydrogenSulfide", pyromat="H2S"),
    "SulfurDioxide": SpeciesRecord("SulfurDioxide", coolprop="SulfurDioxide", pyromat="SO2"),
    "SulfurHexafluoride": SpeciesRecord("SulfurHexafluoride", coolprop="SulfurHexafluoride", pyromat="SF6"),
    "Neon": SpeciesRecord("Neon", coolprop="Neon", pyromat="Ne"),
    "Krypton": SpeciesRecord("Krypton", coolprop="Krypton", pyromat="Kr"),
    "Xenon": SpeciesRecord("Xenon", coolprop="Xenon", pyromat="Xe"),

    # CoolProp-only / real-fluid entries
    "Acetone": SpeciesRecord("Acetone", coolprop="Acetone", pyromat=None),
    "CarbonylSulfide": SpeciesRecord("CarbonylSulfide", coolprop="CarbonylSulfide", pyromat=None),
    "CycloHexane": SpeciesRecord("CycloHexane", coolprop="CycloHexane", pyromat=None),
    "CycloPropane": SpeciesRecord("CycloPropane", coolprop="CycloPropane", pyromat=None),
    "Cyclopentane": SpeciesRecord("Cyclopentane", coolprop="Cyclopentane", pyromat=None),
    "DimethylEther": SpeciesRecord("DimethylEther", coolprop="DimethylEther", pyromat=None),
    "DimethylCarbonate": SpeciesRecord("DimethylCarbonate", coolprop="DimethylCarbonate", pyromat=None),
    "Fluorine": SpeciesRecord("Fluorine", coolprop="Fluorine", pyromat=None),
    "Isohexane": SpeciesRecord("Isohexane", coolprop="Isohexane", pyromat=None),
    "Isopentane": SpeciesRecord("Isopentane", coolprop="Isopentane", pyromat=None),
    "Neopentane": SpeciesRecord("Neopentane", coolprop="Neopentane", pyromat=None),
    "Novec649": SpeciesRecord("Novec649", coolprop="Novec649", pyromat=None),
    "n-Dodecane": SpeciesRecord("n-Dodecane", coolprop="n-Dodecane", pyromat=None),
    "n-Decane": SpeciesRecord("n-Decane", coolprop="n-Decane", pyromat=None),
    "n-Heptane": SpeciesRecord("n-Heptane", coolprop="n-Heptane", pyromat=None),
    "n-Hexane": SpeciesRecord("n-Hexane", coolprop="n-Hexane", pyromat=None),
    "n-Nonane": SpeciesRecord("n-Nonane", coolprop="n-Nonane", pyromat=None),
    "n-Octane": SpeciesRecord("n-Octane", coolprop="n-Octane", pyromat=None),
    "n-Pentane": SpeciesRecord("n-Pentane", coolprop="n-Pentane", pyromat=None),
    "n-Propane": SpeciesRecord("n-Propane", coolprop="n-Propane", pyromat="C3H8"),
})


ALIASES: dict[str, str] = {
    "air": "Air",
    "ar": "Argon",
    "argon": "Argon",
    "co2": "CarbonDioxide",
    "carbon dioxide": "CarbonDioxide",
    "co": "CarbonMonoxide",
    "carbon monoxide": "CarbonMonoxide",
    "he": "Helium",
    "helium": "Helium",
    "h2": "Hydrogen",
    "hydrogen": "Hydrogen",
    "ch4": "Methane",
    "methane": "Methane",
    "n2": "Nitrogen",
    "gn2": "Nitrogen",
    "nitrogen": "Nitrogen",
    "o2": "Oxygen",
    "go2": "Oxygen",
    "lox": "Oxygen",
    "oxygen": "Oxygen",
    "h2o": "Water",
    "steam": "Water",
    "water": "Water",
    "rp-1": "n-Dodecane",
    "rp1": "n-Dodecane",
    "jeta": "n-Dodecane",
    "jet-a": "n-Dodecane",
    "kerosene": "n-Dodecane",
    "octane": "n-Octane",
}

class classproperty(property):

    def __get__(self, obj, owner):
        return self.fget(owner)


class FluidRegistry:

    @staticmethod
    def normalize_name(value: str) -> str:
        return (
            value.strip()
            .lower()
            .replace(" ", "")
            .replace("_", "")
            .replace("-", "")
        )

    @classmethod
    def name(cls, value: str) -> str:
        lookup = cls.normalize_name(value)

        for alias, name in ALIASES.items():
            if cls.normalize_name(alias) == lookup:
                return name

        for name in SPECIES_DATABASE:
            if cls.normalize_name(name) == lookup:
                return name

        raise ValueError(f"Unknown fluid/species name: {value!r}")

    @classmethod
    def record(cls, value: str) -> SpeciesRecord:
        return SPECIES_DATABASE[cls.name(value)]

    @classmethod
    def coolprop_name(cls, value: str) -> str:
        record = cls.record(value)

        if record.coolprop is None:
            raise ValueError(
                f"{record.name!r} is not supported by CoolProp."
            )

        return record.coolprop

    @classmethod
    def pyromat_name(
        cls,
        value: str,
        include_prefix: bool = False,
    ) -> str:
        record = cls.record(value)

        if record.pyromat is None:
            raise ValueError(
                f"{record.name!r} is not supported by PYroMat."
            )

        if include_prefix:
            return f"ig.{record.pyromat}"

        return record.pyromat

    @classmethod
    def supports_coolprop(cls, value: str) -> bool:
        try:
            return cls.record(value).coolprop is not None
        except ValueError:
            return False

    @classmethod
    def supports_pyromat(cls, value: str) -> bool:
        try:
            return cls.record(value).pyromat is not None
        except ValueError:
            return False

    @classmethod
    def supports_both(cls, value: str) -> bool:
        return (
            cls.supports_coolprop(value)
            and cls.supports_pyromat(value)
        )

    @classmethod
    def add_alias(cls, alias: str, name: str) -> None:
        ALIASES[alias] = cls.name(name)

    @classmethod
    def remove_alias(cls, alias: str) -> None:
        ALIASES.pop(alias, None)

    @classproperty
    def names(cls) -> list[str]:
        return sorted(SPECIES_DATABASE.keys())

    @classproperty
    def coolprop_supported_names(cls) -> list[str]:
        return sorted(
            name
            for name, record in SPECIES_DATABASE.items()
            if record.coolprop is not None
        )

    @classproperty
    def pyromat_supported_names(cls) -> list[str]:
        return sorted(
            name
            for name, record in SPECIES_DATABASE.items()
            if record.pyromat is not None
        )

    @classproperty
    def supports_both_names(cls) -> list[str]:
        return sorted(
            name
            for name, record in SPECIES_DATABASE.items()
            if (
                record.coolprop is not None
                and record.pyromat is not None
            )
        )

    @classproperty
    def aliases(cls) -> dict[str, str]:
        return dict(sorted(ALIASES.items()))

    @classmethod
    def show_species(cls) -> list[str]:
        for name in cls.names:
            print(name)

        return cls.names