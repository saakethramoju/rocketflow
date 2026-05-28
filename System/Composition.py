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
        total = sum(state.value for state in self.fraction.values())

        if not np.isclose(total, 1.0, rtol=0.0, atol=1e-6):
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
        # Return the current numeric mass fractions.
        return {
            species: state.value
            for species, state in self.fraction.items()
        }

    @property
    def is_assigned(self) -> bool:
        # Return True if the composition contains at least one species.
        return len(self.fraction) > 0

    def constrain_species(
        self,
        species: str | None = None,
    ) -> None:
        # Select one regular State to adjust so all mass fractions sum to 1.0.
        if not self.is_assigned:
            raise ValueError("Cannot constrain an empty Composition.")

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


    def __getitem__(self, species: str) -> State:
        # Return the species fraction State, or a fixed zero State if absent.
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