from __future__ import annotations

from typing import TYPE_CHECKING

from System import Component

if TYPE_CHECKING:
    from System import Network, State


class Conduction(Component):
    """
    One-dimensional conduction heat transfer between two temperature nodes.

    Positive heat_rate indicates heat transfer from temperature1 to
    temperature2 according to Fourier's law.
    """
    def __init__(self, 
                 name: str, 
                 network: Network,
                 temperature1: State,
                 temperature2: State,
                 thermal_conductivity: State,
                 length: float,
                 conductive_area: float,
                 heat_rate: State | None = None):
        self.setup()

    def evaluate_states(self):
        k = self.thermal_conductivity.value
        A = self.conductive_area.value
        L = self.length.value
        T1 = self.temperature1.value
        T2 = self.temperature2.value

        self.heat_rate.value = k * A / L * (T1 - T2)





class Radiation(Component):
    """
    Diffuse-gray surface radiation exchange between two temperature nodes
    if emissivity is constant at all wavelengths (temperatures).

    Supports arbitrary surface emissivities, radiating areas, and view
    factors. Positive heat_rate indicates net radiative heat transfer
    from temperature1 to temperature2.
    """

    SIGMA = 5.670374419e-8  # W/m^2-K^4

    def __init__(
        self,
        name: str,
        network: Network,
        temperature1: State,
        temperature2: State,
        emissivity1: float,
        emissivity2: float,
        radiative_area1: float,
        radiative_area2: float | None = None,
        view_factor12: float = 1.0,
        heat_rate: State | None = None,
    ):
        self.setup()

        if radiative_area2 is None:
            self.radiative_area2.value = self.radiative_area1.value

    def evaluate_states(self):
        T1 = self.temperature1.value
        T2 = self.temperature2.value

        eps1 = self.emissivity1.value
        eps2 = self.emissivity2.value

        A1 = self.radiative_area1.value
        A2 = self.radiative_area2.value

        F12 = self.view_factor12.value

        if eps1 <= 0.0 or eps1 > 1.0:
            raise ValueError(
                f"{self.name}: emissivity1 must be in (0, 1]. Got {eps1}."
            )

        if eps2 <= 0.0 or eps2 > 1.0:
            raise ValueError(
                f"{self.name}: emissivity2 must be in (0, 1]. Got {eps2}."
            )

        if A1 <= 0.0:
            raise ValueError(
                f"{self.name}: radiative_area1 must be greater than zero. Got {A1}."
            )

        if A2 <= 0.0:
            raise ValueError(
                f"{self.name}: radiative_area2 must be greater than zero. Got {A2}."
            )

        if F12 <= 0.0 or F12 > 1.0:
            raise ValueError(
                f"{self.name}: view_factor12 must be in (0, 1]. Got {F12}."
            )

        denominator = (
            (1.0 - eps1) / (eps1 * A1)
            + 1.0 / (A1 * F12)
            + (1.0 - eps2) / (eps2 * A2)
        )

        if denominator <= 0.0:
            raise ValueError(
                f"{self.name}: invalid radiation denominator ({denominator})."
            )

        self.heat_rate.value = (
            self.SIGMA
            * (T1**4 - T2**4)
            / denominator
        )





class AmbientRadiation(Component):
    """
    Radiation exchange between a surface and a surrounding ambient enclosure.

    Uses the enclosure radiation model with configurable ambient emissivity.
    Positive heat_rate indicates net radiative heat transfer from the solid
    surface to the ambient surroundings.
    """

    SIGMA = 5.670374419e-8  # W/m^2-K^4

    def __init__(
        self,
        name: str,
        network: Network,
        solid_temperature: State,
        ambient_temperature: State | float,
        emissivity: State | float,
        radiative_area: State | float,
        ambient_emissivity: State | float = 1.0,
        heat_rate: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        Ts = self.solid_temperature.value
        Tamb = self.ambient_temperature.value

        eps_s = self.emissivity.value
        eps_amb = self.ambient_emissivity.value

        A = self.radiative_area.value

        if eps_s <= 0.0 or eps_s > 1.0:
            raise ValueError(
                f"{self.name}: emissivity must be in (0, 1]. Got {eps_s}."
            )

        if eps_amb <= 0.0 or eps_amb > 1.0:
            raise ValueError(
                f"{self.name}: ambient_emissivity must be in (0, 1]. Got {eps_amb}."
            )

        if A <= 0.0:
            raise ValueError(
                f"{self.name}: radiative_area must be greater than zero. Got {A}."
            )

        denominator = (
            1.0 / eps_s
            + 1.0 / eps_amb
            - 1.0
        )

        if denominator <= 0.0:
            raise ValueError(
                f"{self.name}: invalid radiation denominator ({denominator})."
            )

        self.heat_rate.value = (
            self.SIGMA
            * A
            * (Ts**4 - Tamb**4)
            / denominator
        )





class Convection(Component):
    """
    Convective heat transfer between a surface and a fluid or ambient node.

    Positive heat_rate indicates heat transfer from surface_temperature to
    fluid_temperature according to Newton's law of cooling.
    """

    def __init__(
        self,
        name: str,
        network: Network,
        surface_temperature: State,
        fluid_temperature: State | float,
        heat_transfer_coefficient: State | float,
        convective_area: State | float,
        heat_rate: State | None = None,
    ):
        self.setup()

    def evaluate_states(self):
        Ts = self.surface_temperature.value
        Tf = self.fluid_temperature.value

        h = self.heat_transfer_coefficient.value
        A = self.convective_area.value

        if h < 0.0:
            raise ValueError(
                f"{self.name}: heat_transfer_coefficient must be nonnegative. Got {h}."
            )

        if A <= 0.0:
            raise ValueError(
                f"{self.name}: convective_area must be greater than zero. Got {A}."
            )

        self.heat_rate.value = h * A * (Ts - Tf)