import numpy as np

from Utilities import FluidRegistry
from .State import State


class Composition:
    """
    Container for species mass-fraction States.
    """

    def __init__(
        self,
        fluid: dict[str, State | float] | str,
    ):
        # Convert a pure species string into a single-species composition.
        if isinstance(fluid, str):
            fluid = {fluid: 1.0}

        # Store all species fractions as mutable State objects.
        self.fraction: dict[str, State] = {
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

    def __getitem__(self, species: str) -> State:
        # Allow composition["O2"] style access.
        return self.fraction[FluidRegistry.name(species)]

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