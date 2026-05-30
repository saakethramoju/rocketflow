import numpy as np

from Utilities import FluidRegistry
from .State import State


class Composition:
    """
    Container for species mass-fraction States.
    """

    def __init__(
        self,
        fluid: dict[str, State | float] | str | None = None,
    ):
        # Initialize empty placeholder composition.
        self.fraction: dict[str, State] = {}
        self._constrained_species: str | None = None
        self._zero_fraction_states: dict[str, State] = {}

        # Empty composition behaves like an unassigned placeholder.
        if fluid is None:
            return

        # Convert a pure species string into a single-species composition.
        if isinstance(fluid, str):
            fluid = {fluid: 1.0}

        # Store all species fractions as mutable State objects.
        self.fraction = {
            FluidRegistry.name(species): (
                value
                if isinstance(value, State)
                else State(float(value))
            )
            for species, value in fluid.items()
        }

        # Require initial mass fractions to sum to 1.0.
        self.validate()


    def validate(
        self,
        atol: float = 1e-6,
    ) -> None:
        # Raise if mass fractions do not sum to 1.
        total = sum(state.value for state in self.fraction.values())

        if not np.isclose(total, 1.0, rtol=0.0, atol=atol):
            raise ValueError(
                f"Composition mass fractions must sum to 1.0. Got {total}."
            )
        

    @property
    def species(self) -> tuple[str, ...]:
        # Return the species names.
        return tuple(self.fraction.keys())

    @property
    def states(self) -> list[State]:
        # Return the mutable fraction States.
        return list(self.fraction.values())

    @property
    def values(self) -> dict[str, float]:
        # Return current numeric mass fractions
        #self.update()
        return {
            species: state.value
            for species, state in self.fraction.items()
        }

    @property
    def is_assigned(self) -> bool:
        # Return True if the composition contains at least one species.
        return len(self.fraction) > 0
        
    def update(self) -> None:
        # Select a constrained species if needed, then enforce sum(Y)=1.
        if self._constrained_species is None:
            self.constrain_species()

        self.enforce_constraint()


    def constrain_species(
        self,
        species: str | None = None,
    ) -> None:
        # Select one regular State to adjust so all mass fractions sum to 1.0.
        if not self.is_assigned:
            return

        if self._constrained_species is not None and species is None:
            return

        if species is None:
            species = next(reversed(self.fraction))

        species = FluidRegistry.name(species)

        if species not in self.fraction:
            raise ValueError(
                f"{species!r} is not present in the composition."
            )

        self._constrained_species = species
        self.enforce_constraint()

    def enforce_constraint(self) -> None:
        # Re-adjust the constrained species if one has been selected.
        if self._constrained_species is None:
            return

        species = self._constrained_species

        value = 1.0 - sum(
            state.value
            for other_species, state in self.fraction.items()
            if other_species != species
        )

        self.fraction[species].value = value


    def copy_from(
        self,
        other: "Composition",
    ) -> None:
        # Copy matching species values from another composition.
        for species in self.species:
            self[species].value = other[species].value


    def __getitem__(self, species: str) -> State:
        # Return a species State, or a fixed zero State if absent.
        #self.update()

        species = FluidRegistry.name(species)

        if species in self.fraction:
            return self.fraction[species]

        if species not in self._zero_fraction_states:
            self._zero_fraction_states[species] = State(0.0)

        return self._zero_fraction_states[species]


    def __iter__(self):
        # Iterate over (species, State) pairs.
        return iter(self.fraction.items())

    def __contains__(self, species: str) -> bool:
        # Check whether a species is present.
        return FluidRegistry.name(species) in self.fraction

    def __len__(self) -> int:
        # Return the number of species.
        return len(self.fraction)
            
    def __or__(self, other: "Composition | None") -> set[str]:
        # Return all unique species from two compositions.
        if other is None:
            return set(self.species)

        return set(self.species) | set(other.species)

    def __and__(
        self,
        other: "Composition | None",
    ) -> tuple[str, ...]:
        # Return shared, unconstrained species between two assigned compositions.
        if (
            other is None
            or not self.is_assigned
            or not other.is_assigned
        ):
            return ()

        return tuple(
            species
            for species in self.species
            if (
                species in other
                and species != self._constrained_species
            )
        )


    def __eq__(self, other) -> bool:
        # Return True if both compositions contain the same species.
        if not hasattr(other, "species"):
            return False

        return set(self.species) == set(other.species)


    def __le__(
        self,
        other: "Composition",
    ) -> bool:
        # Return True if all species in self are present in other.
        return (self - other) == set()
    


    def _new_unchecked(
        self,
        values: dict[str, float],
    ) -> "Composition":
        # Create a Composition result without requiring fractions to sum to 1.
        result = Composition()
        result.fraction = {
            FluidRegistry.name(species): State(float(value))
            for species, value in values.items()
        }
        return result


    def __mul__(
        self,
        scalar: float | int,
    ) -> "Composition":
        # Scale every species fraction by a scalar.
        return self._new_unchecked({
            species: state.value * float(scalar)
            for species, state in self.fraction.items()
        })


    def __rmul__(
        self,
        scalar: float | int,
    ) -> "Composition":
        # Allow scalar * composition.
        return self * scalar


    def __truediv__(
        self,
        scalar: float | int,
    ) -> "Composition":
        # Divide every species fraction by a scalar.
        scalar = float(scalar)

        if scalar == 0.0:
            raise ZeroDivisionError("Cannot divide Composition by zero.")

        return self._new_unchecked({
            species: state.value / scalar
            for species, state in self.fraction.items()
        })


    def __add__(
        self,
        other: "Composition",
    ) -> "Composition":
        # Add like species and keep species that appear in either composition.
        return self._new_unchecked({
            species: self[species].value + other[species].value
            for species in self | other
        })


    def __sub__(
        self,
        other: "Composition",
    ) -> "Composition":
        # Subtract like species; other cannot contain species missing from self.
        missing = set(other.species) - set(self.species)

        if missing:
            raise ValueError(
                "Cannot subtract Composition with species not present in self: "
                f"{sorted(missing)}"
            )

        return self._new_unchecked({
            species: self[species].value - other[species].value
            for species in self.species
        })


    def __str__(self) -> str:
        # Return a compact user-readable composition string.
        if not self.is_assigned:
            return "Composition(<unassigned>)"

        return (
            "Composition("
            + ", ".join(
                f"{species}={state.value:.6g}"
                for species, state in self.fraction.items()
            )
            + ")"
        )

    def __repr__(self) -> str:
        # Return a debug-readable composition string.
        return f"Composition({self.values})"